#
#  This file is part of Sequana software
#
#  Copyright (c) 2016-2026 - Sequana Development Team
#
#  Distributed under the terms of the 3-clause BSD license.
#  The full license is in the LICENSE file, distributed with this software.
#
#  website: https://github.com/sequana/sequana
#  documentation: http://sequana.readthedocs.io
#
##############################################################################
"""Annotation of one sequence for the genome-track viewer (gene track and gene
names in the region table)."""
from urllib.parse import unquote

import numpy as np

__all__ = ["SKIP_TYPES", "annotate_regions", "assign_lanes", "load_annotation", "select_features"]

#: child features and generic containers: the genes and non-coding RNAs are enough for a
#: gene track and this avoids duplicated rows in the region table
SKIP_TYPES = {
    "mRNA",
    "CDS",
    "exon",
    "five_prime_UTR",
    "three_prime_UTR",
    "UTR",
    "start_codon",
    "stop_codon",
    "transcript",
    "region",
    "chromosome",
    "contig",
    "source",
    "supercontig",
    "regulatory",
}

GENBANK_EXTENSIONS = (".gb", ".gbk", ".genbank", ".gbff")


def _gff_attributes(text):
    out = {}
    for item in text.strip().split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            out[key.strip()] = unquote(value.strip())
    return out


def select_features(features):
    """Drop child features and containers (:data:`SKIP_TYPES`). When nothing else is left
    (e.g. a prokaryotic GenBank file made of CDS only), the CDS are the genes."""
    kept = [f for f in features if f["type"] not in SKIP_TYPES]
    return kept or [f for f in features if f["type"] == "CDS"]


def assign_lanes(features):
    """Sort *features* and give each a ``lane`` so that overlapping features are drawn on
    different rows (greedy: first lane where it does not overlap the previous feature).

    :return: the number of lanes used (at least 1)
    """
    features.sort(key=lambda f: (f["start"], f["end"]))
    lane_end = []
    for feature in features:
        for i, end in enumerate(lane_end):
            if feature["start"] > end + 1:
                lane_end[i], feature["lane"] = feature["end"], i
                break
        else:
            feature["lane"] = len(lane_end)
            lane_end.append(feature["end"])
    return max(1, len(lane_end))


def _read_gff3(path, chrom):
    features = []
    names = (chrom, chrom.split(".")[0])
    with open(path) as fin:
        for line in fin:
            if line.startswith("#") or not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9 or cols[0] not in names:
                continue
            attr = _gff_attributes(cols[8])
            ident = attr.get("ID") or attr.get("gene_id") or f"{cols[2]}_{cols[3]}"
            features.append(
                {
                    "start": int(cols[3]),
                    "end": int(cols[4]),
                    "strand": cols[6],
                    "type": cols[2],
                    "name": attr.get("Name") or attr.get("gene_id") or ident,
                    "product": attr.get("product", ""),
                    "id": ident,
                }
            )
    return features


def _read_genbank(path, chrom):
    from sequana.genbank import GenBank

    raw = GenBank(str(path)).genbank_features_parser()
    features = []
    for key in (chrom, chrom.split(".")[0]):
        for item in raw.get(key, []):
            name = item.get("gene") or item.get("gene_id") or item.get("locus_tag") or item["type"]
            features.append(
                {
                    "start": int(item["gene_start"]),
                    "end": int(item["gene_end"]),
                    "strand": item.get("strand") or ".",
                    "type": item["type"],
                    "name": str(name),
                    "product": str(item.get("product") or ""),
                    "id": str(name),
                }
            )
        if features:
            break
    return features


def load_annotation(path, chrom):
    """Features of one sequence from a GFF3 or GenBank file.

    :param path: GFF3 (any extension) or GenBank (``.gb``, ``.gbk``, ``.gbff``) file
    :param chrom: sequence name (the version suffix is ignored, as in sequana_coverage)
    :return: ``(features, lanes)``: dicts with start, end (1-based, inclusive), strand,
        type, name, product, id and lane; and the number of lanes. Child features
        (mRNA, CDS, exon...) are skipped, see :func:`select_features`.
    """
    path = str(path)
    features = _read_genbank(path, chrom) if path.lower().endswith(GENBANK_EXTENSIONS) else _read_gff3(path, chrom)
    features = select_features(features)
    lanes = assign_lanes(features)
    return features, lanes


def annotate_regions(regions, features):
    """For each region, the indices (in *features*) of the features it overlaps.

    :param regions: sequence of dicts with ``start`` and ``end`` (end is exclusive)
    :param features: list of feature dicts, as returned by :func:`load_annotation`
    :return: a list (one entry per region) of lists of indices
    """
    if not features:
        return [[] for _ in regions]
    starts = np.array([f["start"] for f in features])
    ends = np.array([f["end"] for f in features])
    return [[int(i) for i in np.nonzero((starts < r["end"]) & (ends >= r["start"]))[0]] for r in regions]
