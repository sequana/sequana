#  This file is part of Sequana software
#
#  Copyright (c) 2016-2020 - Sequana Development Team
#
#  Distributed under the terms of the 3-clause BSD license.
#  The full license is in the LICENSE file, distributed with this software.
#
#  website: https://github.com/sequana/sequana
#  documentation: http://sequana.readthedocs.io
#
##############################################################################
"""Ribodesigner module"""
import datetime
import json
import math
import re
import shutil
import subprocess
from itertools import chain, product
from pathlib import Path
from typing import Callable, NamedTuple
from urllib.parse import unquote

from sequana import logger, version
from sequana.fasta import FastA
from sequana.lazy import numpy as np
from sequana.lazy import pandas as pd
from sequana.lazy import pylab, pysam
from sequana.melting import tm_probe_rna_dna
from sequana.tools import reverse_complement

logger.setLevel("INFO")

GFF_COLUMNS = ["seqid", "source", "seq_type", "start", "end", "score", "strand", "phase", "attributes"]

# Features that can carry a rRNA annotation, from the most specific (the transcript itself) to the gene. Other types
# (CDS, exon, mRNA, tRNA, ...) are never selected: their product or description often mentions rRNA (e.g. "16S rRNA
# methyltransferase") without being a rRNA.
_RRNA_FEATURE_RANK = {
    "rRNA": 0,
    "rRNA_primary_transcript": 0,
    "misc_RNA": 1,
    "ncRNA": 1,
    "RNA": 1,
    "transcript": 1,
    "rRNA_gene": 2,
    "ncRNA_gene": 2,
    "gene": 2,
}
_RRNA_BIOTYPE_KEYS = ("gene_biotype", "biotype", "transcript_biotype", "ncrna_class")
_RRNA_BIOTYPES = {"rrna", "mt_rrna", "ribosomal_rna"}
_RRNA_TEXT_KEYS = ("product", "description", "note", "name")
# the whole text must describe a rRNA, e.g. "16S ribosomal RNA", "5S rRNA", "large subunit ribosomal RNA", "16S_rRNA",
# "23S ribosomal RNA, partial sequence". "16S rRNA (guanine) methyltransferase" does not match.
_RRNA_TEXT = re.compile(
    r"(?:(?:5|5\.8|16|18|23|25|28)s|ssu|lsu|small subunit|large subunit)?\s*(?:ribosomal rna|rrna)"
    r"(?:\s*\(?\s*(?:5|5\.8|16|18|23|25|28)s\s*\)?)?(?:\s*,?\s*\(?partial(?: sequence)?\)?)?"
)


def parse_gff_attributes(text):
    """Return the attributes of a GFF3 (``key=value;...``) or GTF (``key "value";...``) line as a dictionary.

    Keys are lower-cased and values are URL-decoded.

    >>> parse_gff_attributes("ID=gene:rrsH;biotype=rRNA;description=16S%20ribosomal%20RNA")["description"]
    '16S ribosomal RNA'
    """
    attributes = {}
    if not isinstance(text, str):
        return attributes
    for item in text.split(";"):
        item = item.strip()
        if not item:
            continue
        if "=" in item:
            key, _, value = item.partition("=")
        elif " " in item:
            key, _, value = item.partition(" ")
        else:
            continue
        attributes[key.strip().lower()] = unquote(value.strip().strip('"'))
    return attributes


def is_rrna_attributes(text):
    """Tell whether GFF attributes describe a rRNA, from a biotype or from the product/description/name.

    >>> is_rrna_attributes("ID=gene-rrn5S;gene=rrn5S;gene_biotype=rRNA")
    True
    >>> is_rrna_attributes("ID=rna-x;product=16S ribosomal RNA")
    True
    >>> is_rrna_attributes("ID=gene:b0051;biotype=protein_coding;description=16S rRNA dimethyltransferase")
    False
    """
    attributes = parse_gff_attributes(text)
    if any(attributes.get(key, "").lower() in _RRNA_BIOTYPES for key in _RRNA_BIOTYPE_KEYS):
        return True
    for key in _RRNA_TEXT_KEYS:
        value = attributes.get(key)
        if value and _RRNA_TEXT.fullmatch(value.replace("_", " ").strip().lower()):
            return True
    return False


def find_rrna_like_features(gff):
    """Select the rRNA features of a GFF table whose type is not (only) ``rRNA``.

    Use it when the annotation stores rRNA as ``gene`` or ``ncRNA_gene`` (with a ``gene_biotype``/``biotype`` of
    rRNA or a product/description such as "16S ribosomal RNA"), as ``misc_RNA``, ``ncRNA`` or ``transcript``.
    Only RNA-like feature types are considered (see ``_RRNA_FEATURE_RANK``), so that CDS or genes whose product
    mentions rRNA (e.g. rRNA methyltransferases) are not selected.

    When a gene and its transcript (or several annotations of the same region) are found, only the most specific
    one is kept, so that a rRNA is not counted twice.

    :param gff: DataFrame with the 9 GFF columns (see ``GFF_COLUMNS``).
    :return: the selected rows, sorted by position.
    """
    candidates = gff[gff.seq_type.isin(_RRNA_FEATURE_RANK)]
    candidates = candidates[candidates.attributes.map(is_rrna_attributes)].copy()
    if candidates.empty:
        return candidates
    candidates["_rank"] = candidates.seq_type.map(_RRNA_FEATURE_RANK)
    accepted = []
    for row in candidates.sort_values(["_rank", "seqid", "start"], kind="stable").itertuples():
        overlaps = any(
            a.seqid == row.seqid and a.strand == row.strand and a.start <= row.end and row.start <= a.end
            for a in accepted
        )
        if not overlaps:
            accepted.append(row)
    keep = candidates.loc[[row.Index for row in accepted]].drop(columns="_rank")
    return keep.sort_values(["seqid", "start"])


def _grid(max_probe_len, max_gap):
    """Probe lengths (from ``max_probe_len`` to 41) and gaps (from ``max_gap`` to 11), longest first."""
    return product(range(max_probe_len, 40, -1), range(max_gap, 10, -1))


def _original_candidates():
    return _grid(60, 20)


def _greedy_candidates():
    # 4% of the lengths have no exact fit in 60 x 20. Extending to 70 x 30 and then 80 x 35 gives no failure for
    # sequences up to 200,000 bases.
    return chain(_grid(60, 20), _grid(70, 30), _grid(80, 35))


def _simple_candidates():
    return iter([(50, 15)])


class TilingMethod(NamedTuple):
    """How the probe length and the gap between probes are chosen for a sequence."""

    candidates: Callable  # returns the (probe_len, gap) pairs to try, in order
    exact_fit: bool  # True: the first and last bases must be covered, so the pair must satisfy L = N l + (N-1) g
    short_sequences: bool  # True: sequences shorter than 100 bases get one or two probes (no scan)
    probes_mode: str  # how the - strand probes are placed, see RiboDesigner._get_probes_df


#: Available methods. To add one, write a function yielding (probe_len, gap) pairs and register it here.
TILING_METHODS = {
    "original": TilingMethod(_original_candidates, exact_fit=True, short_sequences=False, probes_mode="generic"),
    "greedy": TilingMethod(_greedy_candidates, exact_fit=True, short_sequences=True, probes_mode="generic"),
    "simple": TilingMethod(_simple_candidates, exact_fit=False, short_sequences=True, probes_mode="simple"),
}


def find_tiling(seq_len, method="greedy"):
    """Probe length and gap to tile a sequence of ``seq_len`` bases.

    For an exact fit, the first probe starts on the first base and the last one ends on the last base:
    ``L = N * probe_len + (N - 1) * gap`` with ``N`` an integer. The first pair of the method that fits is returned.
    Sequences below 50 bases get one 50 nt probe, and sequences below 100 bases get 40 nt probes spaced by 10 nt
    (one or two probes), whatever the method (except ``original``).

    :param seq_len: length of the sequence.
    :param method: one of :data:`TILING_METHODS`.
    :return: (probe_len, gap)
    :raises ValueError: if the method is unknown or no pair fits.

    >>> find_tiling(1542)
    (60, 18)
    >>> find_tiling(116, "greedy")
    (52, 12)
    """
    try:
        tiling = TILING_METHODS[method]
    except KeyError:
        raise ValueError(f"Unknown method '{method}'. Choose from {sorted(TILING_METHODS)}.") from None

    if tiling.short_sequences:
        if seq_len < 50:
            return 50, 15
        elif seq_len < 100:
            return 40, 10

    for probe_len, gap in tiling.candidates():
        if not tiling.exact_fit or (seq_len + gap) % (probe_len + gap) == 0:
            return probe_len, gap

    raise ValueError(f"No probe length / gap combination fits a sequence of {seq_len} bases with method '{method}'.")


class RiboDesigner(object):
    """Design probes for ribosomes depletion.

    From a complete genome assembly FASTA file and a GFF annotation file:

    - Extract genomic sequences corresponding to the selected ``seq_type``.
    - For these selected sequences, design probes computing probe length and inter probe space according to the length of the ribosomale sequence.
    - Detect the highest cd-hit-est identity threshold where the number of probes is inferior or equal to ``max_n_probes``.
    - Report the list of probes in BED and CSV files.

    In the CSV, the oligo names are in column 1 and the oligo sequences in column 2.

    :param fasta: The FASTA file with complete genome assembly to extract ribosome sequences from.
    :param gff: GFF annotation file of the genome assembly. If none provided, assuming the input FastA is
        already made of rRNA.
    :param output_directory: The path to the output directory defaults to ribodesigner.
    :param seq_type: string describing sequence annotation type (column 3 in GFF) to select rRNA from.
    :param max_n_probes: Max number of probes to design
    :param force:  If the `output_directory` already exists, overwrite it.
    :param threads: Number of threads to use in cd-hit clustering.
    :param float identity_step: step to scan the sequence identity (between 0 and 1) defaults to 0.01.
    :param force_clustering:
    :param offtarget_fasta: FASTA file (e.g. the whole genome or transcriptome) used to search for probe
        off-target hits with blastn. Defaults to ``fasta`` when a ``gff`` is provided (hits overlapping the
        selected rRNA features are then ignored). With no GFF, off-target search is only possible when a
        separate ``offtarget_fasta`` is given.
    """

    def __init__(
        self,
        fasta,
        gff=None,
        output_directory="ribodesigner",
        seq_type="rRNA",
        max_n_probes=384,
        force=False,
        threads=4,
        identity_step=0.01,
        force_clustering=False,
        offtarget_fasta=None,
        **kwargs,
    ):
        # Input
        self.fasta = fasta
        self.gff = gff
        self.seq_type = seq_type
        self.max_n_probes = max_n_probes
        self.threads = threads
        self.outdir = Path(output_directory)
        self.identity_step = identity_step
        self.force_clustering = force_clustering
        self.offtarget_fasta = offtarget_fasta if offtarget_fasta else (fasta if gff else None)

        if force:
            self.outdir.mkdir(exist_ok=True)
        else:
            try:
                self.outdir.mkdir()
            except FileExistsError as err:
                logger.error(f"Output directory {output_directory} exists. Use --force or set force=True")
                raise

        # Outputs
        self.filtered_gff = self.outdir / "ribosome_filtered.gff"
        self.ribo_sequences_fasta = self.outdir / "ribosome_sequences.fas"
        self.probes_fasta = self.outdir / "probes_sequences.fas"

        self.clustered_probes_fasta = self.outdir / "clustered_probes.fas"
        self.clustered_probes_csv = self.outdir / "clustered_probes.csv"
        self.clustered_probes_bed = self.outdir / "clustered_probes.bed"
        self.clustered_probes_xlsx = self.outdir / "clustered_probes.xlsx"

        self.unclustered_probes_bed = self.outdir / "unclustered_probes.bed"

        self.output_json = self.outdir / "ribodesigner.json"
        self.probes_report_csv = self.outdir / "probes_report.csv"
        self.offtarget_hits_csv = self.outdir / "offtarget_hits.csv"
        self.coverage_csv = self.outdir / "coverage.csv"

        self.json = {
            "max_n_probes": max_n_probes,
            "identity_step": identity_step,
            "feature": seq_type,
            "method": "greedy",
            "sequana_version": version,
        }

    def get_rna_pos_from_gff(self):
        """Convert a GFF file into a pandas DataFrame filtered according to the
        self.seq_type.
        """
        total_length = 0

        gff = pd.read_csv(self.gff, sep="\t", comment="#", names=GFF_COLUMNS)

        selection = "type"
        if self.seq_type == "auto":
            filtered_gff, selection = find_rrna_like_features(gff), "auto"
        else:
            filtered_gff = gff.query("seq_type == @self.seq_type")
            if filtered_gff.empty and self.seq_type == "rRNA":
                filtered_gff = find_rrna_like_features(gff)
                if len(filtered_gff):
                    selection = "auto"
                    logger.warning(
                        f"No 'rRNA' feature in the annotation: using {len(filtered_gff)} feature(s) annotated as rRNA "
                        f"({', '.join(sorted(filtered_gff.seq_type.unique()))}) from their biotype or product."
                    )
        self.json["feature_selection"] = selection
        self.json["feature_types_used"] = filtered_gff.seq_type.value_counts().to_dict()

        with pysam.Fastafile(self.fasta) as fas:
            with open(self.ribo_sequences_fasta, "w") as fas_out:
                for row in filtered_gff.itertuples():
                    region = f"{row.seqid}:{row.start}-{row.end}"
                    sequence = fas.fetch(region=region)
                    # features on the minus strand must be read as transcribed (RNA orientation)
                    if row.strand == "-":
                        sequence = reverse_complement(sequence)
                    fas_out.write(f">{region}\n{sequence}\n")
                    total_length += len(sequence)
        self.json["input_total_length"] = total_length
        self.json["input_number_sequences"] = filtered_gff.shape[0]

        seq_types = gff.seq_type.unique().tolist()

        self.json["seq_types"] = ",".join(seq_types)
        logger.info(f"Genetic types found in gff: {','.join(seq_types)}")
        logger.info(
            f"Found {filtered_gff.shape[0]} '{self.seq_type}' entries in the annotation file ({total_length}bp long)."
        )
        logger.debug(f"\t" + filtered_gff.to_string().replace("\n", "\n\t"))

        filtered_gff.to_csv(self.filtered_gff)

    def _get_probe_and_step_len(self, seq, method="greedy"):
        """Calculates the probe_len and inter_probe_space for a ribosomal sequence.

        ribo_len = probe_len * n + (inter_probe_space * (n - 1))
        <=>
        n = (ribo_len + inter_probe_space) / (prob_len + inter_probe_space)

        See :func:`find_tiling` and :data:`TILING_METHODS` for the methods.
        """
        try:
            return find_tiling(len(seq.sequence), method)
        except ValueError:
            raise ValueError(
                f"No correct probe length/inter probe space combination was found for {seq.name}"
            ) from None

    @staticmethod
    def _offset_in_reference(seq, ref, k=20, max_offset=100):
        """Position of the first base of ``seq`` in ``ref``, or None if ``seq`` is not a copy of ``ref``.

        The offset is the one supported by most of the k-mers shared by the two sequences (robust to substitutions);
        it must be supported by at least half of the k-mers of ``seq`` and be at most ``max_offset`` nt.
        """
        if abs(len(seq) - len(ref)) > max_offset or len(seq) < 2 * k:
            return None
        index = {}
        for i in range(len(ref) - k + 1):
            index.setdefault(ref[i : i + k], []).append(i)
        votes = {}
        for j in range(len(seq) - k + 1):
            for i in index.get(seq[j : j + k], ()):
                votes[i - j] = votes.get(i - j, 0) + 1
        if not votes:
            return None
        offset, n = max(votes.items(), key=lambda x: (x[1], -abs(x[0])))
        if n < 0.5 * (len(seq) - k + 1) or abs(offset) > max_offset:
            return None
        return offset

    @staticmethod
    def _offset_map(seq, ref, dominant, k=20, window=200, min_run=10):
        """Offset of ``seq`` in ``ref`` along ``ref``, as [(position in ref, offset)], to follow small indels.

        Only the k-mers unique in ``ref`` and shared with ``seq`` on a diagonal within ``window`` nt of the dominant
        offset are used. Runs of fewer than ``min_run`` k-mers (noise) are ignored.
        """
        index = {}
        for i in range(len(ref) - k + 1):
            index.setdefault(ref[i : i + k], []).append(i)
        matches = []
        for j in range(len(seq) - k + 1):
            hits = index.get(seq[j : j + k], ())
            if len(hits) == 1 and abs(hits[0] - j - dominant) <= window:
                matches.append((hits[0], hits[0] - j))
        matches.sort()
        runs = []
        for pos, offset in matches:
            if runs and runs[-1][1] == offset:
                runs[-1][2] += 1
            else:
                runs.append([pos, offset, 1])
        runs = [r for r in runs if r[2] >= min_run]
        points = []
        for pos, offset, _ in runs:
            if not points or points[-1][1] != offset:
                points.append((pos, offset))
        if not points:
            return [(0, dominant)]
        points[0] = (0, points[0][1])
        return points

    @staticmethod
    def _offset_at(offset, pos):
        """Offset at position ``pos`` of the reference; ``offset`` is an int or an offset map."""
        if isinstance(offset, int):
            return offset
        value = offset[0][1]
        for start, off in offset:
            if start > pos:
                break
            value = off
        return value

    def _tiling_references(self, seqs):
        """Place the copies of a same rRNA (e.g. the 16S of each operon) on a common tiling grid.

        The tiling starts at the first annotated base, so copies annotated with a boundary differing by 1 nt would
        get probes shifted by 1 nt, that clustering cannot merge. The longest copy of each group is the reference:
        it sets the probe length and the spacing, and the other copies are tiled on the same grid, aligned on their
        shared k-mers. Identical probes are then removed as duplicates.

        :return: {name: (reference sequence, position of the first base of the sequence in the reference)}
        """
        references = {}
        refs = []
        for seq in sorted(seqs, key=lambda x: -len(x.sequence)):
            for ref in refs:
                offset = self._offset_in_reference(seq.sequence, ref.sequence)
                if offset is not None:
                    offset_map = self._offset_map(seq.sequence, ref.sequence, offset)
                    references[seq.name] = (ref, offset_map[0][1] if len(offset_map) == 1 else offset_map)
                    break
            else:
                refs.append(seq)
                references[seq.name] = (seq, 0)
        shifted = [
            {
                "name": seq.name,
                "reference": ref.name,
                "offset": offset if isinstance(offset, int) else offset[0][1],
                "n_indels": 0 if isinstance(offset, int) else len(offset) - 1,
                "length_difference": len(seq.sequence) - len(ref.sequence),
            }
            for seq in seqs
            for ref, offset in [references[seq.name]]
            if ref is not seq and (offset or len(seq.sequence) != len(ref.sequence))
        ]
        self.json["copies"] = {
            "n_sequences": len(seqs),
            "n_groups": len(refs),
            "n_shifted": len(shifted),
            "shifted": shifted,
        }
        logger.info(f"{len(seqs)} sequences in {len(refs)} group(s) of copies")
        if shifted:
            n_indels = sum(x["n_indels"] for x in shifted)
            logger.warning(
                f"{len(shifted)} copy(ies) of a rRNA differ from the longest copy of their group by their annotated "
                f"boundaries (start shifted by {min(abs(x['offset']) for x in shifted)} to "
                f"{max(abs(x['offset']) for x in shifted)} nt, length difference up to "
                f"{max(abs(x['length_difference']) for x in shifted)} nt) or by {n_indels} indel(s). Their probes were "
                "placed on the grid of that copy so identical probes can be merged. Check the annotation if this is "
                "unexpected."
            )
        return references

    @staticmethod
    def _fit_starts(ref_starts, offset_of, probe_len, length, max_clip=10, anchor_ends=False):
        """Probe starts of the reference grid, moved into a copy; ``offset_of(x)`` gives the shift (nt) to apply to the grid position ``x``.

        A probe sticking out of the copy by at most ``max_clip`` nt is moved inside it so the ends stay covered; the
        others are dropped. With ``anchor_ends``, end-to-end probes are added at an end of the copy that the grid
        leaves uncovered (the copy extends beyond the reference).
        """
        starts = []
        if length < probe_len:
            return starts
        for x in ref_starts:
            x -= offset_of(x)
            if x < -max_clip or x + probe_len > length + max_clip:
                continue
            x = min(max(x, 0), length - probe_len)
            if x not in starts:
                starts.append(x)
        if anchor_ends and starts:
            # cover what the grid leaves at the ends of the copy (an overhang of the copy beyond the reference)
            while min(starts) > 0:
                starts.insert(0, max(0, min(starts) - probe_len))
            while max(starts) < length - probe_len:
                starts.append(min(length - probe_len, max(starts) + probe_len))
        return starts

    def _get_probes_df(self, seq, probe_len, step_len, mode="generic", offset=0, ref_len=None):
        """Generate the Dataframe with probes information.

        Design probes to have end-to-end coverage on the + strand and fill the inter_probe_space present on the + strand with probes designed on the - strand.

        :param seq: A pysam sequence object.
        :param prob_len: The length of the probes calculated by self._get_probe_and_step_len.
        :param step_len: The length of the inter-probe space calculated by self._get_probe_and_step_len.
        :param strand: The strand on which probes are designed.
        :param offset: position of the first base of ``seq`` in the reference copy the tiling is laid on (an int, or
            an offset map when the copy has indels), see
            :meth:`_tiling_references`. The probes of the reference grid that do not fit in ``seq`` are dropped.
        :param ref_len: length of the reference copy (default: the length of ``seq``).
        """
        length = len(seq.sequence)
        ref_len = length if ref_len is None else ref_len
        step = probe_len + step_len

        # + strand probes, on the grid of the reference copy
        ref_starts = list(range(0, ref_len - probe_len + 1, step))
        # sequences below 100 bases get probes that do not reach their end (see find_tiling): add one that does
        if mode != "simple" and ref_starts and ref_starts[-1] + probe_len < ref_len:
            ref_starts.append(ref_len - probe_len)
        is_copy = not (isinstance(offset, int) and offset == 0 and ref_len == length)
        starts = self._fit_starts(
            ref_starts,
            lambda x: self._offset_at(offset, x),
            probe_len,
            length,
            anchor_ends=is_copy and mode != "simple",
        )
        stops = [start + probe_len for start in starts]

        df = pd.DataFrame(
            {
                "name": seq.name,
                "start": starts,
                "stop": stops,
                "strand": "+",
                "score": 0,
            }
        )
        df["sequence"] = [seq.sequence[row.start : row.stop] for row in df.itertuples()]
        df["seq_id"] = df["name"] + f"_+_" + df["start"].astype(str) + "_" + df["stop"].astype(str)

        # - strand probes
        sequence = reverse_complement(seq.sequence)
        # Starts reverse probes to be centered on inter_probe_space of the forward probes
        if mode == "simple":
            ref_rev_starts = ref_starts
        else:
            ref_rev_starts = [int((ref_starts[i + 1] + ref_starts[i]) / 2) for i in range(0, len(ref_starts) - 1)]
        # the reverse complement is numbered from the other end of the copy
        # (position of the probe in the reference as transcribed: ref_len - x - probe_len)
        rev_starts = self._fit_starts(
            ref_rev_starts,
            lambda x: ref_len - length - self._offset_at(offset, ref_len - x - probe_len),
            probe_len,
            length,
        )
        rev_stops = [start + probe_len for start in rev_starts]

        df_rev = pd.DataFrame(
            {
                "name": seq.name,
                "start": rev_starts,
                "stop": rev_stops,
                "strand": "-",
                "score": 0,
            }
        )
        df_rev["sequence"] = [sequence[row.start : row.stop] for row in df_rev.itertuples()]
        df_rev["seq_id"] = df_rev["name"] + f"_-_" + df_rev["start"].astype(str) + "_" + df_rev["stop"].astype(str)

        # Transform to bed coordinates for the reverse_complement
        df_rev["start"] = len(sequence) - df_rev["start"]
        df_rev["stop"] = len(sequence) - df_rev["stop"]
        df_rev.rename(columns={"start": "stop", "stop": "start"}, inplace=True)

        return pd.concat([df, df_rev])

    def get_all_probes(self, method="greedy"):
        """Run all probe design and concatenate results in a single DataFrame."""

        self.json["method"] = method
        probes_dfs = []

        if method not in TILING_METHODS:
            raise ValueError(f"Unknown method '{method}'. Choose from {sorted(TILING_METHODS)}.")
        probes_mode = TILING_METHODS[method].probes_mode

        with pysam.FastxFile(self.ribo_sequences_fasta) as fas:
            seqs = list(fas)
        references = self._tiling_references(seqs)
        for seq in seqs:
            ref, offset = references[seq.name]
            probe_len, step_len = self._get_probe_and_step_len(ref, method)
            probes_dfs.append(
                self._get_probes_df(
                    seq, probe_len, step_len, mode=probes_mode, offset=offset, ref_len=len(ref.sequence)
                )
            )

        self.probes_df = pd.concat(probes_dfs)

        logger.info(f"Number of probes: {len(self.probes_df)}.")
        self.json["n_probes_with_duplicates"] = len(self.probes_df)

        # Probes with the same sequence (the copies of a rRNA operon) are kept in the table, since each has its own
        # position, but only the first is kept for the pool; the others are merged into it at 100% identity.
        df = self.probes_df.reset_index(drop=True)
        duplicated = df.sequence.duplicated()
        first = df.drop_duplicates("sequence").set_index("sequence").seq_id
        df["kept_after_clustering"] = ~duplicated
        df["cluster_representative"] = df.sequence.map(first)
        df["cluster_identity"] = 100.0
        df["bed_color"] = df.kept_after_clustering.map({True: "21,128,0", False: "128,64,0"})
        self.probes_df = df
        logger.info(f"Number of unique probes: {int((~duplicated).sum())}.")

    def export_all_probes_to_fasta(self):
        """From the self.probes_df, export to FASTA and CSV files."""

        with open(self.probes_fasta, "w") as fas:
            for row in self.probes_df.itertuples():
                fas.write(f">{row.seq_id}\n{row.sequence}\n")

    def _n_unique_probes(self):
        """Number of probes with a distinct sequence."""
        return int(self.probes_df.sequence.nunique())

    def clustering_needed(self, force=False):
        """Checks if a clustering is needed.

        :param force: force clustering even if unecessary.
        """

        # Do not cluster if number of probes already inferior to defined threshold
        if force or self.force_clustering:
            logger.warning(f"You specified the --force-clustering option. Performing clustering")
            return True
        elif not (force or self.force_clustering) and self._n_unique_probes() <= self.max_n_probes:
            logger.info(
                f"Number of probes {self._n_unique_probes()} already inferior to {self.max_n_probes}. No clustering will be performed."
            )
            return False
        else:
            return True

    @staticmethod
    def _cdhit_word_size(identity):
        """Word size (-n) recommended by the cd-hit-est documentation for an identity threshold."""
        if identity >= 0.9:
            return 8
        elif identity >= 0.88:
            return 7
        elif identity >= 0.85:
            return 6
        elif identity >= 0.8:
            return 5
        return 4

    @staticmethod
    def _parse_clstr(path):
        """Parse a cd-hit ``.clstr`` file into {member: (representative, identity_percent)}."""
        members, cluster = {}, []
        rep_name = None

        def flush():
            for name, pct in cluster:
                members[name] = (rep_name, pct)

        with open(path) as fin:
            for line in fin:
                if line.startswith(">Cluster"):
                    if cluster:
                        flush()
                    cluster, rep_name = [], None
                    continue
                match = re.search(r">(.+?)\.\.\. (\*|at [+-]/([\d.]+)%)", line)
                if not match:
                    continue
                name = match.group(1)
                if match.group(2) == "*":
                    rep_name = name
                    cluster.append((name, 100.0))
                else:
                    cluster.append((name, float(match.group(3))))
        if cluster:
            flush()
        return members

    def cluster_probes(self):
        """Use cd-hit-est to cluster highly similar probes."""
        logger.info("Clustering probes")
        outdir = (
            Path(self.clustered_probes_fasta).parent
            / f"cd-hit-est-{datetime.datetime.today().isoformat(timespec='seconds', sep='_')}"
        )
        outdir.mkdir()
        log_file = outdir / "cd-hit.log"

        res_dict = {"seq_id_thres": [], "n_probes": []}

        for seq_id_thres in np.arange(0.8, 1, self.identity_step).round(3):
            tmp_fas = outdir / f"clustered_{seq_id_thres}.fas"
            cmd = (
                f"cd-hit-est -i {self.probes_fasta} -o {tmp_fas} -c {seq_id_thres} "
                f"-n {self._cdhit_word_size(seq_id_thres)} -T {self.threads} -d 0"
            )
            logger.debug(f"Clustering probes with command: {cmd} (log in '{log_file}').")

            with open(log_file, "a") as f:
                subprocess.run(cmd, shell=True, check=True, stdout=f)

            res_dict["seq_id_thres"].append(seq_id_thres)
            res_dict["n_probes"].append(len(FastA(tmp_fas)))

        # Add number of probes without clustering
        res_dict["seq_id_thres"].append(1)
        res_dict["n_probes"].append(self._n_unique_probes())

        df = pd.DataFrame(res_dict)
        # Extract the best identity threshold
        # If we force the clustering, the number of probes is already below max_n_probes even without clustering
        # (seq_id_thres==1), so we need to take the min of the max_b_probes and the current number of probes
        # with no clustering
        M = min(self.max_n_probes, self._n_unique_probes() - 1)
        best_thres = df.query("n_probes <= @M").seq_id_thres.max()

        # n_probes corresponding to the best identity threshold
        # best_thres is NaN when even the lowest identity does not bring the number of probes under the budget
        n_probes = (
            df.query("seq_id_thres == @best_thres").loc[:, "n_probes"].values[0]
            if not np.isnan(best_thres)
            else self._n_unique_probes()
        )
        self.json["n_probes"] = int(n_probes)
        self.json["best_thres"] = best_thres
        self.json["results"] = res_dict

        # Dataframe with number of probes for each cdhit identity threshold
        df = pd.DataFrame(self.json["results"])

        if not np.isnan(best_thres):
            logger.info(f"Best clustering threshold: {best_thres}, with {n_probes} probes.")
            shutil.copy(outdir / f"clustered_{best_thres}.fas", self.clustered_probes_fasta)
            kept_probes = [seq.name for seq in FastA(outdir / f"clustered_{best_thres}.fas")]
            self.probes_df["kept_after_clustering"] = self.probes_df.seq_id.isin(kept_probes)
            self.probes_df["bed_color"] = self.probes_df.kept_after_clustering.map(
                {True: "21,128,0", False: "128,64,0"}
            )
            self.probes_df["clustering_thres"] = best_thres
            members = self._parse_clstr(outdir / f"clustered_{best_thres}.fas.clstr")
            self.probes_df["cluster_representative"] = self.probes_df.seq_id.map(
                lambda x: members.get(x, (x, 100.0))[0]
            )
            self.probes_df["cluster_identity"] = self.probes_df.seq_id.map(lambda x: members.get(x, (x, 100.0))[1])
        else:
            logger.warning(
                f"No identity threshold was found to have as few as {self.max_n_probes} probes. Keep all probes. Set a valid value with --max-n-probes between {df.n_probes.min()} (min) and {df.n_probes.max()} (max)"
            )
            with open(self.clustered_probes_fasta, "w") as fas:
                for row in self.probes_df.query("kept_after_clustering == True").itertuples():
                    fas.write(f">{row.seq_id}\n{row.sequence}\n")

        self.clustering_df = df.sort_values("seq_id_thres")

        return self.probes_df.query("kept_after_clustering == True")

    def plot(self):

        # Dataframe with number of probes for each cdhit identity threshold
        pylab.clf()
        df = pd.DataFrame(self.json["results"])

        import seaborn as sns  # local import to speed up imports

        p = sns.lineplot(data=df, x="seq_id_thres", y="n_probes", markers=["o"])
        p.axhline(
            self.max_n_probes,
            alpha=0.8,
            linestyle="--",
            color="red",
            label="max number of probes requested",
        )
        pylab.xlabel("Sequence identity", fontsize=16)
        pylab.ylabel("Number of probes", fontsize=16)

        best_thres = self.json["best_thres"]
        n_probes = self.json["n_probes"]

        if not np.isnan(best_thres):
            pylab.plot(best_thres, n_probes, "o", label="Final number of probes")
            pylab.legend()
        else:
            logger.warning(
                f"No identity threshold was found to have as few as {self.max_n_probes} probes. Keep all probes. Set a valid value with --max-n-probes between {df.n_probes.min()} (min) and {df.n_probes.max()} (max)"
            )
            with open(self.clustered_probes_fasta, "w") as fas:
                for row in self.probes_df.query("kept_after_clustering == True").itertuples():
                    fas.write(f">{row.seq_id}\n{row.sequence}\n")

    def export_to_csv_bed(self):
        """Export final results to CSV and BED files"""

        # first we save unclustered probes as a BED file
        self.probes_df.to_csv(
            self.unclustered_probes_bed,
            sep="\t",
            index=False,
            header=None,
            columns=[
                "name",
                "start",
                "stop",
                "sequence",
                "score",
                "strand",
                "start",
                "stop",
                "bed_color",
            ],
        )

        # then we save final CSV and BED file.
        if self.clustering_needed():
            df = self.cluster_probes()
        else:
            df = self.probes_df.query("kept_after_clustering == True")
            with open(self.clustered_probes_fasta, "w") as fas:
                for row in df.itertuples():
                    fas.write(f">{row.seq_id}\n{row.sequence}\n")

        df.to_csv(self.clustered_probes_csv, index=False, columns=["seq_id", "sequence"])

        # xlsx for IDT format
        dd = df[["seq_id", "sequence"]].copy()
        dd.columns = ["Pool name", "Sequence"]
        dd["Pool name"] = "poolOne"
        dd.to_excel(self.clustered_probes_xlsx, index=False)

        # save clustered probes in BED format
        df.to_csv(
            self.clustered_probes_bed,
            sep="\t",
            index=False,
            header=None,
            columns=[
                "name",
                "start",
                "stop",
                "sequence",
                "score",
                "strand",
                "start",
                "stop",
                "bed_color",
            ],
        )

    # ------------------------------------------------------------------ probe QC
    @staticmethod
    def _gc_content(sequence):
        """GC fraction (0-1) of a sequence."""
        sequence = sequence.upper()
        return (sequence.count("G") + sequence.count("C")) / len(sequence) if sequence else 0.0

    @staticmethod
    def _tm(sequence, na_molar=0.05):
        """Approximate melting temperature (Celsius), salt-adjusted GC formula for long oligos.

        Tm = 81.5 + 16.6 log10([Na+]) + 41 (G+C)/N - 675/N

        Fallback for probes with ambiguous bases (see :meth:`_tm_nn`). This is an estimate (no nearest-neighbour
        model, no RNA:DNA correction): use it to compare probes within a pool, not as an absolute value.
        """
        n = len(sequence)
        if n == 0:
            return float("nan")
        gc = RiboDesigner._gc_content(sequence)
        return 81.5 + 16.6 * math.log10(na_molar) + 41.0 * gc - 675.0 / n

    @staticmethod
    def _tm_nn(sequence, na_mm=100, probe_nm=1000):
        """Melting temperature (Celsius) of the DNA probe bound to its RNA target (nearest-neighbour model).

        Uses the RNA/DNA hybrid parameters of Sugimoto et al. (1995), see :func:`sequana.melting.tm_probe_rna_dna`.
        Conditions: ``na_mm`` mM Na+ and ``probe_nm`` nM probe in large excess over the target. Returns None if the
        sequence contains ambiguous bases (e.g. N), for which the model is not defined.
        """
        try:
            return tm_probe_rna_dna(sequence, na_mm=na_mm, probe_nm=probe_nm)
        except ValueError:
            return None

    @staticmethod
    def _longest_stem(sequence, min_loop=None):
        """Longest stretch whose reverse complement is also present in the sequence.

        With ``min_loop=None`` the two copies may overlap (what matters for a self-dimer). With an integer, the
        two copies must be separated by at least ``min_loop`` bases (what matters for a hairpin).
        """
        seq = sequence.upper()
        n = len(seq)
        for k in range(n // 2 if min_loop is not None else n, 0, -1):
            for i in range(n - k + 1):
                rc = reverse_complement(seq[i : i + k])
                start = seq.find(rc, i + k + min_loop) if min_loop is not None else seq.find(rc)
                if start != -1:
                    return k
        return 0

    @staticmethod
    def _longest_homopolymer(sequence):
        """Length of the longest run of a single nucleotide."""
        runs = [len(m.group(0)) for m in re.finditer(r"(.)\1*", sequence.upper())]
        return max(runs) if runs else 0

    def compute_qc(self, gc_range=(0.3, 0.7), max_homopolymer=5, max_hairpin=8, max_selfdimer=10):
        """Add quality-control columns to ``self.probes_df``.

        Columns added: ``length``, ``gc``, ``tm``, ``max_homopolymer``, ``hairpin_stem``, ``selfdimer_stem`` and
        ``qc_flags`` (comma-separated list of ``low_gc``, ``high_gc``, ``homopolymer``, ``hairpin``,
        ``selfdimer``; empty when the probe passes).

        ``tm`` is the RNA:DNA nearest-neighbour Tm (see :meth:`_tm_nn`). For a probe with ambiguous bases, for
        which the model is not defined, a salt-adjusted GC estimate is used instead (see :meth:`_tm`); the number
        of such probes is stored in ``json["qc"]["n_tm_estimated"]``.
        ``hairpin_stem`` and ``selfdimer_stem`` are the longest self-complementary stretches (hairpins need a loop
        of at least 3 bases); they are sequence-only heuristics, not free-energy predictions.

        :param gc_range: (min, max) acceptable GC fraction.
        :param max_homopolymer: longest acceptable single-nucleotide run.
        :param max_hairpin: longest acceptable hairpin stem (bases).
        :param max_selfdimer: longest acceptable self-dimer stretch (bases).
        """
        df = self.probes_df
        df["length"] = df.sequence.str.len()
        df["gc"] = df.sequence.map(self._gc_content).round(3)

        tm_nn = df.sequence.map(self._tm_nn)
        df["tm"] = tm_nn.fillna(df.sequence.map(self._tm)).round(1)
        tm_method = "nearest-neighbour RNA/DNA (Sugimoto 1995), 100 mM Na+, 1 uM probe"

        df["max_homopolymer"] = df.sequence.map(self._longest_homopolymer)
        df["hairpin_stem"] = df.sequence.map(lambda x: self._longest_stem(x, min_loop=3))
        df["selfdimer_stem"] = df.sequence.map(self._longest_stem)

        def flags(row):
            out = []
            if row.gc < gc_range[0]:
                out.append("low_gc")
            if row.gc > gc_range[1]:
                out.append("high_gc")
            if row.max_homopolymer > max_homopolymer:
                out.append("homopolymer")
            if row.hairpin_stem > max_hairpin:
                out.append("hairpin")
            if row.selfdimer_stem > max_selfdimer:
                out.append("selfdimer")
            return ",".join(out)

        df["qc_flags"] = df.apply(flags, axis=1)

        self.json["qc"] = {
            "gc_range": list(gc_range),
            "max_homopolymer": max_homopolymer,
            "max_hairpin": max_hairpin,
            "max_selfdimer": max_selfdimer,
            "tm_method": tm_method,
            "n_tm_estimated": int(tm_nn.isna().sum()),
            "n_flagged": int((df.qc_flags != "").sum()),
            "gc_mean": float(df.gc.mean()),
            "tm_mean": float(df.tm.mean()),
            "tm_min": float(df.tm.min()),
            "tm_max": float(df.tm.max()),
        }
        logger.info(
            f"QC: {self.json['qc']['n_flagged']} / {len(df)} probes flagged ({gc_range=}, {max_homopolymer=}, {max_hairpin=}, {max_selfdimer=})."
        )

    # ------------------------------------------------------------ off-target
    def _rrna_intervals(self):
        """Return {seqid: [(start, end), ...]} (1-based, inclusive) of the selected rRNA features."""
        intervals = {}
        if self.gff:
            gff = pd.read_csv(self.filtered_gff, index_col=0)
            for row in gff.itertuples():
                intervals.setdefault(row.seqid, []).append((int(row.start), int(row.end)))
        return intervals

    def check_offtarget(self, min_identity=80, min_aligned=25, evalue=10):
        """Search probes with blastn against ``offtarget_fasta`` and report hits outside the rRNA loci.

        Hits on both strands are considered: probes are designed on both strands, and the strand of the
        target is not known in advance, so this is conservative. Hits overlapping a selected rRNA feature are
        on-target and discarded.

        Adds ``n_offtarget`` (number of off-target hits) and ``best_offtarget_identity`` to
        ``self.probes_df``, and writes all hits to ``offtarget_hits.csv``.

        :param min_identity: minimum percent identity of a reported hit.
        :param min_aligned: minimum alignment length (bases) of a reported hit.
        :param evalue: blastn e-value cutoff (permissive by default because probes are short).
        """
        if not self.offtarget_fasta:
            logger.warning("No off-target FASTA available (no GFF and no offtarget_fasta): skipping off-target check.")
            return
        if shutil.which("blastn") is None or shutil.which("makeblastdb") is None:
            raise RuntimeError("blastn/makeblastdb not found in PATH; required for the off-target check.")

        db = self.outdir / "offtarget_db" / "db"
        db.parent.mkdir(exist_ok=True)
        log_file = self.outdir / "blast.log"
        with open(log_file, "a") as log:
            subprocess.run(
                ["makeblastdb", "-in", str(self.offtarget_fasta), "-dbtype", "nucl", "-out", str(db)],
                check=True,
                stdout=log,
                stderr=log,
            )
            cols = "qseqid sseqid pident length qstart qend sstart send evalue bitscore"
            res = subprocess.run(
                [
                    "blastn",
                    "-task",
                    "blastn-short",
                    "-query",
                    str(self.probes_fasta),
                    "-db",
                    str(db),
                    "-perc_identity",
                    str(min_identity),
                    "-evalue",
                    str(evalue),
                    "-num_threads",
                    str(self.threads),
                    "-outfmt",
                    f"6 {cols}",
                ],
                check=True,
                capture_output=True,
                text=True,
            )
            log.write(res.stderr)

        hits = pd.DataFrame([line.split("\t") for line in res.stdout.splitlines()], columns=cols.split()).astype(
            {
                "pident": float,
                "length": int,
                "sstart": int,
                "send": int,
                "evalue": float,
                "bitscore": float,
            }
        )
        hits = hits.query("length >= @min_aligned")

        intervals = self._rrna_intervals()

        def on_target(row):
            lo, hi = sorted((row.sstart, row.send))
            return any(lo <= e and hi >= s for s, e in intervals.get(row.sseqid, []))

        if len(hits):
            hits = hits[~hits.apply(on_target, axis=1)]
        hits.to_csv(self.offtarget_hits_csv, index=False)

        counts = hits.groupby("qseqid").size()
        best = hits.groupby("qseqid").pident.max()
        self.probes_df["n_offtarget"] = self.probes_df.seq_id.map(counts).fillna(0).astype(int)
        self.probes_df["best_offtarget_identity"] = self.probes_df.seq_id.map(best)

        self.json["offtarget"] = {
            "min_identity": min_identity,
            "min_aligned": min_aligned,
            "n_probes_with_offtarget": int((self.probes_df.n_offtarget > 0).sum()),
            "n_hits": int(len(hits)),
        }
        logger.info(
            f"Off-target: {self.json['offtarget']['n_probes_with_offtarget']} / {len(self.probes_df)} probes "
            f"with at least one hit outside rRNA loci."
        )

    # -------------------------------------------------------------- coverage
    def compute_coverage(self):
        """Report how much of each rRNA sequence is covered by the probes.

        Clustering merges similar probes (for instance the copies of a rRNA operon), so a probe removed at one
        position is still represented by a kept probe at another position. Two coverages are therefore reported
        for each sequence:

        - ``covered_designed_pct``: positions spanned by at least one designed probe (before clustering), on both
          strands together (``covered_plus_pct`` / ``covered_minus_pct`` for each strand alone). This is the
          coverage of the final pool, since every removed probe is represented by a kept probe at least
          ``min_merged_identity`` percent identical.
        - ``covered_direct_pct``: positions spanned by a probe kept at that position.

        Gaps (stretches with no designed probe) are listed per sequence. Results go to ``coverage.csv`` and
        ``json["coverage"]``. Probe coordinates are relative to the rRNA sequence as transcribed (see
        :meth:`get_rna_pos_from_gff`).
        """
        lengths = {}
        with pysam.FastxFile(str(self.ribo_sequences_fasta)) as fas:
            for seq in fas:
                lengths[seq.name] = len(seq.sequence)

        def pct(sub, name, length):
            cov = np.zeros(length, dtype=int)
            for row in sub[sub.name == name].itertuples():
                cov[max(0, int(row.start)) : min(length, int(row.stop))] += 1
            return cov

        def gaps(cov):
            out, start = [], None
            for i, c in enumerate(cov):
                if c == 0 and start is None:
                    start = i
                elif c > 0 and start is not None:
                    out.append([start, i])
                    start = None
            if start is not None:
                out.append([start, len(cov)])
            return out

        df = self.probes_df
        kept = df[df.kept_after_clustering]
        records = []
        for name, length in lengths.items():
            cov = pct(df, name, length)
            g = gaps(cov)
            records.append(
                {
                    "name": name,
                    "length": length,
                    "covered_designed_pct": round(100 * float((cov > 0).mean()), 2),
                    "covered_direct_pct": round(100 * float((pct(kept, name, length) > 0).mean()), 2),
                    "covered_plus_pct": round(100 * float((pct(df[df.strand == "+"], name, length) > 0).mean()), 2),
                    "covered_minus_pct": round(100 * float((pct(df[df.strand == "-"], name, length) > 0).mean()), 2),
                    "n_gaps": len(g),
                    "largest_gap": max((b - a for a, b in g), default=0),
                    "gaps": g,
                }
            )

        total = sum(lengths.values())

        def overall(key):
            return round(sum(r[key] * r["length"] for r in records) / total, 2) if total else 0.0

        merged = df[~df.kept_after_clustering]
        self.coverage_df = pd.DataFrame(records)
        self.coverage_df.assign(gaps=self.coverage_df.gaps.map(json.dumps)).to_csv(self.coverage_csv, index=False)
        self.json["coverage"] = {
            "covered_designed_pct": overall("covered_designed_pct"),
            "covered_direct_pct": overall("covered_direct_pct"),
            "covered_plus_pct": overall("covered_plus_pct"),
            "covered_minus_pct": overall("covered_minus_pct"),
            "n_gaps": int(sum(r["n_gaps"] for r in records)),
            "n_merged_probes": int(len(merged)),
            "min_merged_identity": float(merged.cluster_identity.min()) if len(merged) else 100.0,
            "per_sequence": records,
        }
        logger.info(
            f"Coverage: {self.json['coverage']['covered_designed_pct']}% of the rRNA covered by the pool "
            f"({self.json['coverage']['covered_direct_pct']}% directly by kept probes, the rest through "
            f"probes merged at >= {self.json['coverage']['min_merged_identity']}% identity)."
        )

    def export_probes_report(self):
        """Write the full probes table (coordinates, QC, off-target, clustering status) to CSV."""
        self.probes_df.to_csv(self.probes_report_csv, index=False)

    def export_to_json(self):
        with open(self.output_json, "w") as fout:
            json.dump(self.json, fout, indent=4, sort_keys=True)

    def run(self, method="greedy", offtarget=False):
        """Run the whole design.

        :param method: probe/step length scanning method (original, greedy, simple).
        :param offtarget: if True, search probes for off-target hits (requires blastn, see :meth:`check_offtarget`).
        """
        if self.gff:
            self.get_rna_pos_from_gff()
        else:
            shutil.copy(self.fasta, self.ribo_sequences_fasta)

        self.get_all_probes(method=method)
        self.compute_qc()
        self.export_all_probes_to_fasta()
        if offtarget:
            self.check_offtarget()
        self.export_to_csv_bed()
        self.compute_coverage()
        self.export_probes_report()
        self.export_to_json()
