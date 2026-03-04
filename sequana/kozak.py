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
import functools
import os
import random
import re
from collections import Counter, defaultdict
from contextlib import contextmanager

import colorlog
from tqdm import tqdm

from sequana.fasta import FastA
from sequana.gff3 import GFF3
from sequana.lazy import numpy as np
from sequana.lazy import pandas as pd
from sequana.lazy import pylab
from sequana.tools import reverse_complement

logger = colorlog.getLogger(__name__)


class Motif:
    def __init__(self, motif):
        self.df = motif

    def _get_entropy(self):
        from scipy.stats import entropy

        freq = self.df[["A", "C", "G", "T"]].to_numpy()
        freq = freq / freq.sum(axis=1, keepdims=True)
        IC = 2 - np.array([entropy(row, base=2) for row in freq])
        return IC

    entropy = property(_get_entropy)

    def plot_entropy(self):
        from pylab import axvline, plot, ylabel, ylim

        axvline(0.5, c="k")
        plot(self.df.index, self.entropy)
        ylabel("Bits")
        ylim([0, 1.5])


class Kozak:
    """

        k = KOSAC("fasta", "gff", feature="gene")

    Filter only to keep ATG since others seems to ncRNA

    - raw Kozak sequence names and counts
    - a Kozak is e.g GGCRGG  . first position is the less important

    - for the enumeration of kmers, get of the rid of the Ns

    - odds ratio have 4 cases depending on the on enumeration:
        use entire genome
        use chromosome by chromosome
        use of gene on genome
        use gene on chromosomes


    Table of counts of Kozak sequences without dna ambiguities.
    - across the entire genome
    - by chromosomes

        counts = k.get_all_kmer_counts()
        counts_chroms = k.get_all_kmer_counts_by_chromosome()
        counts_genes = k.get_all_kmer_counts_genes_only()

        # proportions of kmer in genes:
        sum(list(counts_genes.values())) / length_genome

        # counts in chroms should equal counts in genomes:
        Sgenes = sum([sum(list(counts_chroms[x].values())) for x in counts_chroms.keys()])
        Sgenome = sum(list(counts.values()))

    ::

        k = Kozak("ecoli_MG1655.fa", "ecoli_MG1655.gff", "gene", "ID")
        df = k.get_data()
        k.plot_logo(df.query("start_codon=='ATG'"))


    """

    def __init__(self, fasta, gff, genetic_type="gene", attribute="ID", light=True):
        self.fasta = FastA(fasta)
        if isinstance(gff, GFF3):
            self.gff = gff
        else:
            self.gff = GFF3(gff, light=light)
        logger.info("scanning GFF")
        self._valid_genetic_types = self.gff.features  # some overhead but required
        self.genetic_type = genetic_type
        self.attribute = attribute
        self.set_context()

        # internal boundaries to precompute data
        self.RIGHT = 100
        self.LEFT = 100

        # place holder for various metrics
        self.metrics = {}

    @property
    def genetic_type(self):
        return getattr(self, "_genetic_type", None)

    @genetic_type.setter
    def genetic_type(self, value):
        if value not in self._valid_genetic_types:
            raise ValueError(f"Invalid genetic_type '{value}'. " f"Must be one of {sorted(self._valid_genetic_types)}")
        self._genetic_type = value

    def set_context(self, left_kozak=6, right_kozak=6, keep_ATG_only=True, include_start_codon=False):
        assert left_kozak > 0
        assert right_kozak > 0
        self._left_kozak = left_kozak
        self._right_kozak = right_kozak
        self._keep_ATG_only = keep_ATG_only
        self._include_start_codon = include_start_codon
        self.get_random_atg_contexts.cache_clear()
        try:
            del self._cached_df
        except AttributeError:
            pass

    def _get_left_kozak(self):
        return self._left_kozak

    left_kozak = property(_get_left_kozak)

    def _get_right_kozak(self):
        return self._right_kozak

    right_kozak = property(_get_right_kozak)

    def _get_keep_ATG_only(self):
        return self._keep_ATG_only

    def _set_keep_ATG_only(self, value):
        assert value in [True, False]
        self._keep_ATG_only = value

    keep_ATG_only = property(_get_keep_ATG_only, _set_keep_ATG_only)

    def _get_include_start_codon(self):
        return self._include_start_codon

    def _set_include_start_codon(self, value):
        assert value in [True, False]
        self._include_start_codon = value

    include_start_codon = property(_get_include_start_codon, _set_include_start_codon)

    def _compute(self):
        # Store all kozak sequences with generous left and right values
        if hasattr(self, "_cached_df"):
            return self._cached_df

        genetic_type = self.genetic_type
        LEFT = self.LEFT
        RIGHT = self.RIGHT

        logger.info("1. Checking consistency between FastA and GFF files")
        if len(set(self.gff.contig_names).intersection(self.fasta.names)) != len(self.fasta.names):
            logger.warning("GFF and FASTA sequence identifiers have different length")

        # trick to not load all data in memory but only selected genetic type
        self.gff.skip_types = [x for x in self._valid_genetic_types if x != self.genetic_type]
        self.gff._df = None

        # Reading GFF
        gff = self.gff.df

        # we split by chrom to get the sequence one by one.
        data = []
        warning_message = None
        for chrom, subdf in gff.groupby("seqid"):
            sequence = self.fasta[chrom].upper()
            for row in subdf.itertuples(index=False):
                try:
                    ID = getattr(row, self.attribute)
                except AttributeError:
                    ID = ""
                    if warning_message is None:  # show message only once
                        warning_message = (
                            f"Attribute {self.attribute} not found in GFF file. Setting ID to empty string."
                        )
                        logger.warning(warning_message)

                strand = row.strand
                start = row.start
                end = row.stop

                if strand == "+":
                    codon = sequence[start - 1 : start - 1 + 3]
                    kozak_left = sequence[start - 1 - LEFT : start - 1]
                    kozak_right = sequence[start + 3 : start + 3 + RIGHT]
                    data.append([chrom, ID, strand, start, end, codon, kozak_left, kozak_right])
                elif strand == "-":
                    codon = sequence[end - 3 : end]
                    codon = reverse_complement(codon)
                    koz = sequence[end : end + LEFT]
                    kozak_left = reverse_complement(koz)
                    koz = sequence[end - RIGHT - 3 : end - 3]
                    kozak_right = reverse_complement(koz)
                    data.append([chrom, ID, strand, start, end, codon, kozak_left, kozak_right])

        df = pd.DataFrame(data)
        try:
            df.columns = ["chrom", "ID", "strand", "start", "end", "start_codon", "kozak_left", "kozak_right"]
        except ValueError:  # empty dataframe
            pass

        # Cleanup region where selected left and right sub sequence do not
        # have the corret length. This could be a gene right at the border of contig.
        N0 = len(df)
        df = df[[len(x) == LEFT for x in df["kozak_left"].values]]
        df = df[[len(x) == RIGHT for x in df["kozak_right"].values]]
        ratio = len(df) / N0
        self.metrics["feature"] = self.genetic_type
        logger.info(
            f"Filtered to keep only rows with correct left and right Kozak lengths: {len(df)} / {N0} ({ratio:.1%})"
        )
        self.metrics["ATG_ratio"] = sum(df["start_codon"] == "ATG") / len(df)
        self._cached_df = df

        return df

    def get_data(self):
        # do not touch the original
        df = self._compute().copy()

        N = len(df)
        if self.keep_ATG_only:
            df = df[df["start_codon"] == "ATG"]
            logger.info(f"Filtered to keep only ATG start codons: {len(df)} / {N} ({len(df)/N:.1%})")
            self.atg_contribution = 100.0 * len(df) / N
        else:
            n = len(df[df["start_codon"] == "ATG"])
            self.atg_contribution = 100 * float(n) / N
        df["kozak_left"] = [x[-self._left_kozak :] for x in df["kozak_left"]]
        df["kozak_right"] = [x[: self._right_kozak] for x in df["kozak_right"]]
        df["sequence"] = [x + y + z for x, y, z in zip(df["kozak_left"], df["start_codon"], df["kozak_right"])]
        freqs = Counter(df["sequence"])
        df["frequency"] = [freqs[x] / len(df) for x in df["sequence"].values]

        return df

    def filter_dataframe(self, df, strand=None, query=None, genes_set=None, attribute=None):
        if strand == "strand+":
            df = df.query("strand == '+'")
        elif strand == "strand-":
            df = df.query("strand == '-'")

        if query:
            df = df.query(query)

        if genes_set is not None:
            df = df[df[attribute].isin(genes_set)]

        return df

    @functools.lru_cache(maxsize=1)  # None means no limit on cache size
    def get_gc_per_chromosome(self, quiet=True):
        GCs = []
        chrom_names = []
        # we use a list rather than a dictionary to keep same order as in the
        # fasta
        from sequana.tools import fast_gc_content

        for chrom in tqdm(self.fasta.names, disable=quiet):
            sequence = self.fasta.sequences[self.fasta.names.index(chrom)]
            GCs.append(fast_gc_content(sequence))
            chrom_names.append(chrom)
        del sequence
        return chrom_names, GCs

    def plot_GC_per_chromosome(self, ylim=[0, 100]):
        chrom_names, GCs = self.get_gc_per_chromosome()
        pylab.plot(range(0, len(GCs)), [100 * x for x in GCs], "o-")
        if len(GCs) < 100:
            pylab.xticks(range(0, len(GCs)), chrom_names, rotation=90)
        pylab.xlabel("Chromosome")
        pylab.grid()
        pylab.ylim(ylim)
        pylab.ylabel("GC (%)")

    def _get_logo_data(self, df=None):

        if df is None:
            df = self.get_data()

        try:
            Nl = len(df["kozak_left"].iloc[0])
            Nr = len(df["kozak_right"].iloc[0])
        except KeyError:
            N = self.left_kozak
            return {"status": "Warning", "msg": "No data found"}

        if self.include_start_codon:
            pos = [defaultdict(int) for x in range(Nl + 3 + Nr)]
            for seql, start, seqr in zip(df.kozak_left.values, df.start_codon.values, df.kozak_right.values):
                # left sequence
                for i in range(Nl):
                    pos[i][seql[i]] += 1
                # start codon
                for i in range(3):
                    pos[Nl + i][start[i]] += 1

                # right sequence
                for i in range(Nr):
                    pos[Nl + 3 + i][seqr[i]] += 1
        else:
            pos = [defaultdict(int) for x in range(Nl + Nr)]
            for seql, seqr in zip(df.kozak_left.values, df.kozak_right.values):
                for i in range(Nl):
                    pos[i][seql[i]] += 1

                for i in range(Nr):
                    pos[Nl + i][seqr[i]] += 1

        logo_data = pd.DataFrame(pos).fillna(0)
        if "N" in logo_data.columns:
            del logo_data["N"]

        logo_data = logo_data.divide(logo_data.sum(axis=1), axis=0)

        L, R = self.left_kozak, self.right_kozak
        if self.include_start_codon:
            logo_data.index = list(range(-L, 0)) + list(range(1, R + 4))
        else:
            logo_data.index = list(range(-L, 0)) + list(range(4, R + 4))

        return logo_data

    def plot_logo(self, df=None, ax=None, color_scheme="colorblind"):
        assert color_scheme in ["colorblind", "classic"], "color_scheme must be colorblind or classic"
        if color_scheme == "colorblind":
            color_scheme = {
                "A": "#0072B2",
                "C": "#E69F00",
                "G": "#009E73",
                "T": "#D55E00",
            }
        elif color_scheme == "classic":
            color_scheme = {"A": "green", "C": "blue", "T": "red", "G": "orange"}

        logo_data = self._get_logo_data(df)
        self._plot_logo(logo_data, ax=ax, color_scheme=color_scheme)

        return logo_data

    def _add_purine_pyrimidine(self, df):
        df["R"] = df["A"] + df["G"]
        df["Y"] = df["C"] + df["T"]
        return df

    def plot_logo_purine_pyrimidine(self, df=None, ax=None):
        """
        df is the output of :meth:`get_data`
        """
        logo_data = self._get_logo_data(df)
        logo_data = self._add_purine_pyrimidine(logo_data)
        self._plot_logo(logo_data[["R", "Y"]], ax=ax, color_scheme={"Y": "purple", "R": "#78bc00"})
        return logo_data

    def get_entropy(self, motif):
        return Motif(motif).entropy

    @contextmanager
    def temporary_lr(self, left=None, right=None):
        old_left = self.left_kozak
        old_right = self.right_kozak
        try:
            if left is not None:
                self._left_kozak = left
            if right is not None:
                self._right_kozak = right
            yield
        finally:
            self._left_kozak = old_left
            self._right_kozak = old_right

    def _get_KL_data(self, df=None, n_boot=500, ci=95, left=None, right=None):

        if left or right:
            self.get_random_atg_contexts.cache_clear()
            with self.temporary_lr(left=left, right=right):
                logo_data = self._get_logo_data(df=df)
                background = self.get_random_atg_contexts()
                mu, low, high = self.bootstrap(
                    logo_data,
                    background,
                    n_boot=n_boot,
                    ci=ci,
                )
        else:
            logo_data = self._get_logo_data(df=df)
            background = self.get_random_atg_contexts()
            mu, low, high = self.bootstrap(
                logo_data,
                background,
                n_boot=n_boot,
                ci=ci,
            )

        df = pd.DataFrame({"position": logo_data.index, "KL_divergence": mu, "ci_low": low, "ci_high": high})
        return df

    def plot_KL_divergence(self, df=None, ax=None, n_boot=500, ci=95):
        KL = self._get_KL_data(df, n_boot, ci)

        if ax is None:
            fig, ax = pylab.subplots(figsize=(12, 4))
        ax.plot(KL["position"], KL["KL_divergence"], marker="o")
        ax.fill_between(KL["position"], KL["ci_low"], KL["ci_high"], color="gray", alpha=0.3)
        ax.set_xlabel("Position relative to start codon")
        ax.set_ylabel("KL divergence (bits)")
        ax.set_title("Positional KL divergence vs background distribution")
        return KL

    def _plot_logo(self, data, color_scheme=None, ax=None):
        import logomaker

        data = data.copy()
        indices = data.index
        data.reset_index(inplace=True, drop=True)

        logo = logomaker.Logo(data, ax=ax, color_scheme=color_scheme, stack_order="fixed")
        for x in [0.25, 0.5, 0.75]:
            pylab.axhline(x, color="grey", zorder=-1, ls="--")
        pylab.axvline(self._left_kozak - 0.5, color="k", lw=2)

        if self.include_start_codon:
            pylab.axvline(self._left_kozak + 2.5, color="k", lw=2)
        pylab.xticks(range(len(data)), indices)
        _ = pylab.yticks([])

    def plot_kozak_chi2(self, motif=None, GC_mode="context", fontsize=10, noplot=False):
        """
        for GC computation, 3 options:
            1. genomic ATG background excluding annotated starts.go)od for GC bias, codon bias, local sequence structure.
            2. genome wide base composition. simple and fast but inflates signal in GC biases genomes. does not control for ATG specific context. This is a simple a genome-wide background is estimated assuming strand symmetry:
                G = C = GC / 2
                A = T = AT / 2
            3. uniform background. comparable across species, no computation but biologically naive and misleading for GC genomes

        """
        if motif is None:
            motif = self._get_logo_data()
        df = motif  # just an alias

        if GC_mode == "context":
            atg = self.get_random_atg_contexts()
            seq = atg["kozak_left"].str.cat(atg["kozak_right"])
            GC = seq.str.count("[GC]").sum() / seq.str.len().sum()
        elif GC_mode == "genome":
            GC = self.metrics.get("GC", self.fasta.GC_content())
        elif GC_mode == "uniform":
            GC = 0.5
        else:
            raise ValueError("GC_mode must be context, genome, uniform")

        # computes chi2
        from scipy.stats import chisquare

        expected = {"A": 100 * (1 - GC) / 2.0, "C": 100 * GC / 2.0, "G": 100 * GC / 2.0, "T": 100 * (1 - GC) / 2}
        observed = 100 * df[["A", "C", "G", "T"]].copy()

        chi2 = []
        for _, row in observed.iterrows():
            chi2_stat, p_value = chisquare(
                [row["A"], row["C"], row["G"], row["T"]], f_exp=[expected[n] for n in ["A", "C", "G", "T"]]
            )
            chi2.append(p_value)
        observed["chi2"] = chi2
        if noplot:
            return observed

        df = observed
        fig, ax = pylab.subplots(figsize=(12, 10))
        # Plotting values
        ax.matshow(df[["A", "C", "G", "T"]].T, cmap="viridis", vmin=0, vmax=50)

        # Adding text annotations
        df = df[["A", "C", "G", "T"]].T
        for i in range(df.shape[0]):
            for j in range(df.shape[1]):
                ax.text(j, i, f"{df.values[i, j]:.1f}", ha="center", va="center", color="red", fontsize=fontsize)
        pylab.yticks([0, 1, 2, 3], ["A", "C", "G", "T"])
        stars = []
        for p in observed["chi2"]:
            if p < 0.01:
                stars.append("**")
            elif p < 0.05:
                stars.append("*")
            else:
                stars.append(" ")

        pylab.xticks(range(len(observed)), stars, fontsize=fontsize)
        pylab.axvline(self.left_kozak - 0.5, lw=2, color="k")
        if self.include_start_codon:
            pylab.axvline(self._left_kozak + 2.5, color="k", lw=2)

        return observed

    @functools.lru_cache(maxsize=1)
    def _get_annotated_starts_by_chrom(self):
        """
        Returns:
            dict[str, set[int]]: chromosome -> set of start positions (1-based)
        """

        self.gff.skip_types = [x for x in self.gff.features if x != self.genetic_type]
        df = self.gff.df

        starts = defaultdict(set)
        for chrom, start in zip(df.seqid.values, df.start.values):
            starts[chrom].add(int(start))

        return starts

    def _is_annotated_start(self, chromosome, position):
        return position in self._get_annotated_starts_by_chrom().get(chromosome, set())

    @functools.lru_cache(maxsize=None)  # None means no limit on cache size
    def get_random_atg_contexts(self, Nmax=10000, quiet=True):
        # user defined values assuming users know what they are doing

        #
        chrom_lengths = [l for l in self.fasta.lengths]
        cum_lengths = np.cumsum(chrom_lengths)
        total_length = cum_lengths[-1]

        lefts = []
        rights = []
        starts = []

        logger.info("Get annotated start position")
        annotated = self._get_annotated_starts_by_chrom()  # must be precomputed set

        logger.info("Scanning reference")
        while len(lefts) < Nmax:

            # random global position
            gpos = random.randint(0, total_length - 3)

            # map global position to chromosome
            chrom_idx = np.searchsorted(cum_lengths, gpos, side="right")
            if chrom_idx == 0:
                local_pos = gpos
            else:
                local_pos = gpos - cum_lengths[chrom_idx - 1]

            seq = self.fasta.sequences[chrom_idx]
            chrom = self.fasta.names[chrom_idx]

            # 2 jumps to next ATG
            local_pos = seq.find("ATG", local_pos)
            if local_pos == -1:
                local_pos = seq.find("atg", local_pos)
                if local_pos == -1:
                    continue

            # 3 context checks
            if local_pos + 1 in annotated.get(chrom, set()):
                continue

            left = seq[local_pos - self.left_kozak : local_pos]
            right = seq[local_pos + 3 : local_pos + 3 + self.right_kozak]

            if len(left) != self.left_kozak or len(right) != self.right_kozak:
                continue

            lefts.append(left.upper())
            rights.append(right.upper())
            starts.append("ATG")

        return pd.DataFrame({"kozak_left": lefts, "start_codon": starts, "kozak_right": rights})

    def export_meme(self, filename, name="Kozak"):
        """PWM compatible with standard motif scanners"""
        df = self.get_data()
        pwm = self._get_logo_data(df)[["A", "C", "G", "T"]]
        with open(filename, "w") as f:
            f.write("MEME version 4\n\n")
            f.write("ALPHABET= ACGT\n\n")
            f.write("strands: + -\n\n")
            f.write(f"MOTIF {name}\n")
            f.write(f"letter-probability matrix: alength= 4 w= {len(pwm)}\n")
            for _, row in pwm.iterrows():
                f.write(" ".join(f"{row[x]:.6f}" for x in "ACGT") + "\n")

    def kl_vs_random_atg(self, motif_df, random_df):
        """
        Compute position-wise divergence using Kullback–Leibler (KL) divergence
        between Kozak contexts and random (non-annotated) ATG contexts.

        This quantifies how specific the Kozak signal is compared
        to generic ATG neighborhoods.

        Compute positional KL divergence between the observed Kozak motif
        and a background nucleotide distribution.

        This method computes, for each position i:

            D_KL(P_i || Q) = sum_i (P_i x log2(P_i / Q_i))

        where P_i is the observed nucleotide frequency distribution
        at position i, and Q is a fixed background distribution.

        Notes
        -----
        - This quantity is mathematically related to Shannon entropy.
        - When Q is uniform, this is equivalent to classical
          sequence logo information content.
        - The result is deterministic for a given motif_df and  random_df input pair

        Parameters
        ----------
        motif_df : pandas.DataFrame
            DataFrame with columns ['A', 'C', 'G', 'T'] containing
            nucleotide frequencies per position.

        background : background distribution.

        Returns
        -------
        numpy.ndarray
            KL divergence (bits) for each motif position.
        """
        from scipy.stats import entropy

        kl = []
        for (_, p_row), (_, q_row) in zip(
            motif_df[["A", "C", "G", "T"]].iterrows(), random_df[["A", "C", "G", "T"]].iterrows()
        ):
            p = p_row.values.astype(float) + 1e-12
            q = q_row.values.astype(float) + 1e-12
            kl.append(entropy(p, q, base=2))

        return np.asarray(kl)

    def bootstrap(self, df, contexts, n_boot=500, ci=95, sample_size=200):

        # order ACGT
        n = len(contexts)
        if self.include_start_codon:
            seq = contexts["kozak_left"].str.cat(contexts["start_codon"]).str.cat(contexts["kozak_right"])
            L = self.left_kozak + 3 + self.right_kozak
        else:
            seq = contexts["kozak_left"].str.cat(contexts["kozak_right"])
            L = self.left_kozak + self.right_kozak

        boot = np.zeros((n_boot, L))

        for b in tqdm(range(n_boot)):
            idx = np.random.randint(0, n, min(n, sample_size))
            sample = contexts.iloc[idx]
            bg = self._get_logo_data(sample)
            boot[b] = self.kl_vs_random_atg(df, bg)

        mean = boot.mean(axis=0)
        low = np.percentile(boot, (100 - ci) / 2, axis=0)
        high = np.percentile(boot, 100 - (100 - ci) / 2, axis=0)

        return mean, low, high

    def get_KSI(self, KL):
        return sum(KL.query("position<0 and position>=-6")["KL_divergence"]) / 6

        # just 20 boostrap is enough to get an idea
        # KL = self._get_KL_data(left=100, right=6, n_boot=20)
        # and 50bp away by 50bp should give us the noise.
        # bg = KL.iloc[0:50]["KL_divergence"]
        # mu, sigma = mean(bg), std(bg)
        # the signal
        # KL = k._get_KL_data()


class KLAnalysis:
    def __init__(self, df):
        self.KL = df

    def get_KSI(self):
        """average across Kozak length (6bp before ATG)

        Average across 6 bp (Kozak sequence)
        """
        return sum(self.KL.query("position<0 and position>=-6")["KL_divergence"]) / 6

    def get_total_information(self, min=-1e6, max=0):
        """Area under the curve for position<0

        This removes dilution from averaging.
        Independent of window scaling. Measures total constraint.
        """
        return self.KL.query("position<@max and position>=@min")["KL_divergence"].sum()

    def get_peak_strength(self):
        """max peak strength"""
        return self.KL["KL_divergence"].max()

    def get_peak_position(self):
        """max peak position"""
        return self.KL["KL_divergence"].argmax()

    def get_signal_concentration(self):
        # peak concentration measures how concentrated is the signal. sharp pike == high C (eukar), broader = lower C (proka)
        Ipeak = self.get_peak_strength()
        Itot = self.get_total_information()
        signal_concentration = Ipeak / Itot
        return signal_concentration

    def get_III(self):
        # Alternative Single-Value Metric that balances global and local structure:
        from math import sqrt

        Ipeak = self.get_peak_strength()
        Itot = self.get_total_information()
        return sqrt(Ipeak * Itot)

    def get_power(self):
        # power metric: emphasize sharpness
        Ipeak = self.get_peak_strength()
        Itot = self.get_total_information()
        P = Ipeak * Ipeak / Itot
        return P

    def compute_W50(self):
        df_up = self.KL.query("position < 0")
        KL_values = df_up["KL_divergence"].values

        I_total = KL_values.sum()
        if I_total == 0:
            return 0

        sorted_KL = np.sort(KL_values)[::-1]
        cumulative = np.cumsum(sorted_KL)

        W50 = np.searchsorted(cumulative, 0.5 * I_total) + 1
        return W50


class ConsensusBuilder:
    def __init__(self, df):

        self.df = df.copy()

    def _iupac(self, bases):
        from sequana.iuapc import dna_ambiguities_r

        return dna_ambiguities_r.get(bases, "N")

    def _information_content(self, row):
        """Shannon information (DNA max = 2 bits)."""
        p = row.values
        p = p[p > 0]
        H = -np.sum(p * np.log2(p))
        return 2 - H

    def get_consensus(
        self, mode="majority", threshold=0.25, relative=0.8, strong=0.6, majority=0.5, info_threshold=1.0
    ):
        """
        Parameters
        ----------
        mode : str
            majority | threshold | relative | information | max_only
        threshold : float
            used in threshold mode
        relative : float
            keep bases >= relative * max_frequency
        strong : float
            uppercase if max >= strong
        info_threshold : float
            uppercase if information >= this value
        """

        consensus = []

        for _, row in self.df.iterrows():

            if mode == "max_only":
                base = row.idxmax()
                consensus.append(base)
                continue

            elif mode == "majority":
                if row.max() > majority:
                    base = row.idxmax()
                    consensus.append(base.upper())
                else:
                    consensus.append("n")

            elif mode == "threshold":
                bases = row[row >= threshold].index.tolist()
                letter = self._iupac(bases)

                if row.max() >= strong:
                    consensus.append(letter.upper())
                else:
                    consensus.append(letter.lower())

            elif mode == "relative":
                m = row.max()
                bases = row[row >= relative * m].index.tolist()
                letter = self._iupac(bases)

                if m >= strong:
                    consensus.append(letter.upper())
                else:
                    consensus.append(letter.lower())

            elif mode == "information":
                info = self._information_content(row)
                m = row.max()

                # select dominant bases (relative strategy)
                bases = row[row >= 0.8 * m].index.tolist()
                letter = self._iupac(bases)

                if info >= info_threshold:
                    consensus.append(letter.upper())
                else:
                    consensus.append(letter.lower())

            else:
                raise ValueError(f"Unknown mode: {mode}")

        return "".join(consensus)

    def all_consensus(
        self, mode="majority", threshold=0.25, relative=0.8, strong=0.6, majority=0.5, info_threshold=1.0
    ):

        for mode in ["majority", "threshold", "relative", "information", "max_only"]:
            res = self.get_consensus(
                mode=mode,
                threshold=threshold,
                relative=relative,
                strong=strong,
                majority=majority,
                info_threshold=info_threshold,
            )
            res = res
            print(f"{mode}: {res}")


class KozakAddon(Kozak):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    @functools.lru_cache(maxsize=None)  # None means no limit on cache size
    def builddata(self):
        dm, dp = self.get_kmer_counts()

        df = pd.DataFrame({"plus": dp, "minus": dm})
        df = df.sum(axis=1)
        df = df.fillna(0).sort_values(ascending=False)

        self._df = df
        return df

    def _get_df(self):
        if self._df is None:
            self.builddata()
        return self._df

    df = property(_get_df)

    def find_kmers(self, sequence, pattern=r"GCC[AG]CC"):
        """
        >>> k.find_kmers("GCCACC")
        True
        >>> k.find_kmers("AAAAAA")
        False
        """

        kozak_pattern = re.compile(pattern)
        match = kozak_pattern.search(sequence)
        if match:
            return True
        else:
            return False

    @functools.lru_cache(maxsize=None)  # None means no limit on cache size
    def get_all_kmer_counts(self, k=6, reverse=False, quiet=True):
        """Get all kmers from the entire genome"""
        counts = defaultdict(int)
        for sequence in tqdm(self.fasta.sequences, disable=quiet):
            count = self._get_kmer_from_sequence(sequence, k=k, reverse=reverse)

            for x, y in count.items():
                counts[x] += y
        return counts

    def _get_kmer_from_sequence(self, sequence, k=6, reverse=False):
        if reverse:
            sequence = reverse_complement(sequence)
        counts = defaultdict(int)
        for i in range(0, len(sequence) - k + 1):
            seq = sequence[i : i + k]
            if "N" not in seq:
                counts[sequence[i : i + k]] += 1
        return counts

    @functools.lru_cache(maxsize=None)  # None means no limit on cache size
    def get_all_kmer_counts_genes_only(self, k=6, genetic_type="gene", reverse=False):
        counts = defaultdict(int)

        df = self.gff.df.query("genetic_type==@genetic_type and strand=='+'")

        # here we loop through all genes
        sequences = self.fasta.sequences
        names = self.fasta.names

        for index, row in df.iterrows():
            # we find the sequence corresponding to the seqid
            try:
                index_seq = names.index(row.seqid)
                # and the given gene
                sequence = sequences[index_seq][row.start - 1 : row.stop]
                # to extract all kmers inside the gene
                count = self._get_kmer_from_sequence(sequence, k=k, reverse=reverse)

                for x, y in count.items():
                    counts[x] += y
            except ValueError:
                pass
        return counts

    @functools.lru_cache(maxsize=None)  # None means no limit on cache size
    def get_all_kmer_counts_by_chromosome(self, k=6, reverse=False):
        counts = {}
        for chrom, sequence in tqdm(zip(self.fasta.names, self.fasta.sequences)):
            count = self._get_kmer_from_sequence(sequence, k=k, reverse=reverse)
            counts[chrom] = count
        return counts

    """def get_all_kmers_normalised(self, k=6):
        kmers = self.get_all_kmer_counts_by_chromosome(k=k)
        kmers = pd.DataFrame(kmers)
        # kmers = kmers.divide(self.fasta.get_lengths_as_dict())*1000000
        return kmers
    """

    @functools.lru_cache(maxsize=None)  # None means no limit on cache size
    def get_odd_ratio(self, mode="all"):
        odds = []

        counts = self.df
        Ngenes = self.df.sum()
        N_kmers = len(self.df)
        print(f"Number of genes : {Ngenes}")
        print(f"Number of unique kmers : {N_kmers}")

        if mode == "all":
            kmers = self.get_all_kmer_counts()
            genome_size = sum(list(kmers.values()))
        elif mode == "gene":
            kmers = self.get_all_kmer_counts_genes_only()
            genome_size = sum(list(kmers.values()))
        elif mode == "chromosome":
            pass

        odds = []
        for kmer in counts.index:
            count = counts[kmer]

            A = count / Ngenes
            B = kmers[kmer] / genome_size
            odds.append([kmer, count, kmers[kmer], A, B, A / B])

        odds = pd.DataFrame(odds)
        odds.columns = ["kmer", "count", "count_genome", "ratio1", "ratio2", "odds"]
        return odds

    def plot_scatter_odds_ratio_gene_vs_genome(self):
        # odds ratio are computed on the same set of kmer but compared
        # to different random kmer distribution
        odds = self.get_odd_ratio(mode="all")
        odds_gene = self.get_odd_ratio(mode="gene")

        # so we can plot a scatter plot
        # odds['count;[ and odds_gene['count'] are the same
        pylab.scatter(odds["odds"], odds_gene["odds"], c=odds["count"], alpha=0.5, s=4 * odds["count"])
        pylab.xlabel("odds ratio (vs all entire genome)")
        pylab.ylabel("odds ratio (vs all entire coding genes)")
        pylab.colorbar()

    def plot_scatter_odds_ratio_annotated(self, pattern=r"GCC[AG]CC"):
        odds = self.get_odd_ratio(mode="gene")
        found = [10 if self.find_kmers(x, pattern=pattern) else 0 for x in odds["kmer"]]
        odds["found"] = found

        print(odds.query("found==10"))
        pylab.clf()
        pylab.scatter(odds["count"], odds["odds"], c=odds["found"], alpha=0.5, s=4 * odds["odds"])
        pylab.xlabel("Count")
        pylab.ylabel("odds")
        pylab.colorbar()
        return odds

    def _get_kmer_counts_strand(self, strand):
        df = self.get_data()
        chroms = df.chrom.unique()
        kozak = defaultdict(int)

        for chrom in tqdm(chroms):
            for seq in df.query("chrom==@chrom and strand==@strand").kozak_left:
                kozak[seq] += 1

        return kozak

    def get_kmer_counts(self):
        """Get kmer counts for both strands. Returns (minus_strand_counts, plus_strand_counts)"""
        dm = self._get_kmer_counts_strand("-")
        dp = self._get_kmer_counts_strand("+")
        return dm, dp

    def plot_logo_all_kmers(self):
        counts = self.get_all_kmer_counts()

        pos = [defaultdict(int) for x in range(6)]
        for sequence, N in counts.items():
            for i, letter in enumerate(sequence):
                pos[i][letter] += N

        logo_data = pd.DataFrame(pos)
        dd = logo_data.divide(logo_data.sum(axis=1), axis=0)
        self._plot_logo(dd)
        return dd

    def plot_cumulated(self):
        cs = np.cumsum(self.df)
        mid = cs.max() / 2
        pylab.plot(cs.values)
        pylab.axhline(mid)
        pylab.xlabel("Number of unique 6-mers")
        pylab.ylabel("Number of genes")
