#
#  This file is part of Sequana software
#
#  Copyright (c) 2016-2022 - Sequana Development Team
#
#  Distributed under the terms of the 3-clause BSD license.
#  The full license is in the LICENSE file, distributed with this software.
#
#  website: https://github.com/sequana/sequana
#  documentation: http://sequana.readthedocs.io
#
##############################################################################
import functools
from collections import defaultdict

import colorlog
from tqdm import tqdm

from sequana.fasta import FastA
from sequana.gff3 import GFF3
from sequana.lazy import numpy as np
from sequana.lazy import pandas as pd
from sequana.lazy import pylab
from sequana.tools import reverse_complement

logger = colorlog.getLogger(__name__)


class Kozak:
    """

        k = KOSAC("fasta", "gff", feature="gene")


    Filter only to keep ATG since others seems to ncRNA

    - raw Kozak sequence names and counts
    - a Kozak is e.g GGCRGG  . first position is the less important

    - for the enueration of kmers, get of the rid of the Ns

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

        k = Kozak("ecoli_MG1655.fa", "ecoli_MG1655.gff")
        df = k.get_data("CDS", attr="gene")
        k.plot_logo(df.query("start_codon=='ATG'"))


    """

    def __init__(self, fasta, gff, genetic_type="gene", attribute="ID"):
        self.fasta = FastA(fasta)
        self.gff = GFF3(gff)
        self.gff.get_attributes()
        self.genetic_type = genetic_type
        self.attribute = attribute

        assert attribute in self.gff.get_attributes()
        assert genetic_type in self.gff.features

        self._left_kozak = 6
        self._right_kozak = 3
        self._keep_ATG_only = True
        self._df = None

    def _get_left_kozak(self):
        return self._left_kozak

    left_kozak = property(_get_left_kozak)

    def _get_right_kozak(self):
        return self._right_kozak

    right_kozak = property(_get_right_kozak)

    def _get_keep_ATG_only(self):
        return self._keep_ATG_only

    keep_ATG_only = property(_get_keep_ATG_only)

    @functools.lru_cache(maxsize=None)  # None means no limit on cache size
    def builddata(self):
        plus = self.get_strand_plus_data()
        minus = self.get_strand_minus_data()

        df = pd.DataFrame({"plus": plus, "minus": minus})
        # df = df.fillna(0).sort_values(by="plus", ascending=False)
        df = df.sum(axis=1)
        df = df.fillna(0).sort_values(ascending=False)

        self._df = df
        return df

    def _get_df(self):
        if self._df is None:
            self.builddata()
        return self._df

    df = property(_get_df)

    def plot_cumulated(self):
        cs = np.cumsum(self.df)
        mid = cs.max() / 2
        pylab.plot(cs.values)
        pylab.axhline(mid)
        pylab.xlabel("Number of unique 6-mers")
        pylab.ylabel("Number of genes")

    @functools.lru_cache(maxsize=None)  # None means no limit on cache size
    def get_data(self):

        genetic_type = self.genetic_type
        attr = self.attribute

        # we split by chrom to get the sequence one by one.
        data = []
        for chrom, sequence in tqdm(zip(self.fasta.names, self.fasta.sequences)):

            subdf = self.gff.df.query("genetic_type==@genetic_type and seqid==@chrom")

            for index, row in subdf.iterrows():
                ID = row[attr]
                strand = row["strand"]
                start = row["start"]
                end = row["stop"]
                if strand == "+":
                    codon = sequence[start - 1 : start - 1 + 3]
                    if self.keep_ATG_only:
                        if codon == "ATG":
                            kozak_left = sequence[start - 1 - self._left_kozak : start - 1]
                            kozak_right = sequence[start + 3 : start + 3 + self._right_kozak]
                            data.append([chrom, ID, strand, start, end, codon, kozak_left, kozak_right])
                    else:  # we keep every gene
                        kozak_left = sequence[start - 1 - self._left_kozak : start - 1]
                        kozak_right = sequence[start - 1 + 3 : start - 1 + 3 + self._right_kozak]
                        data.append([chrom, ID, strand, start, end, codon, kozak_left, kozak_right])
                elif strand == "-":
                    codon = sequence[end - 3 : end]
                    codon = reverse_complement(codon)
                    if self.keep_ATG_only:
                        if codon == "ATG":
                            koz = sequence[end : end + self._left_kozak]
                            kozak_left = reverse_complement(koz)
                            koz = sequence[end - self._right_kozak - 3 : end - 3]
                            kozak_right = reverse_complement(koz)
                            data.append([chrom, ID, strand, start, end, codon, kozak_left, kozak_right])
                    else:  # we keep all genes
                        koz = sequence[end : end + self._left_kozak]
                        kozak_left = reverse_complement(koz)
                        koz = sequence[end - self.right_kozak - 3 : end - 3]
                        kozak_right = reverse_complement(koz)
                        data.append([chrom, ID, strand, start, end, codon, kozak_left, kozak_right])
        df = pd.DataFrame(data)
        try:
            df.columns = ["chrom", "ID", "strand", "start", "end", "start_codon", "kozak_left", "kozak_right"]
        except ValueError:  # empty dataframe
            pass

        df["sequence"] = [x + y + z for x, y, z in zip(df["kozak_left"], df["start_codon"], df["kozak_right"])]

        from collections import Counter

        freqs = Counter(df["sequence"])
        df["frequency"] = [freqs[x] / len(df) for x in df["sequence"].values]

        return df

    def get_kmer_counts_plus_strand(self):
        """From the full set, extract count on strand + only"""
        return self._get_kmer_counts_strand("+")

    def get_kmer_counts_minus_strand(self):
        """From the full set, extract count on strand - only"""
        return self._get_kmer_counts_strand("-")

    def _get_kmer_counts_strand(self, strand):
        df = self.get_data()
        chroms = df.chrom.unique()
        kozak = defaultdict(int)

        for chrom in tqdm(chroms):
            for seq in df.query("chrom==@chrom and strand==@strand").kozak:
                kozak[seq] += 1

        return kozak

    def get_kmer_counts(self):
        dm = self._get_kmer_counts_strand("-")
        dp = self._get_kmer_counts_strand("-")

    @functools.lru_cache(maxsize=None)  # None means no limit on cache size
    def get_gc_per_chromosome(self):
        GCs = []
        chrom_names = []
        # we use a list rather than a dictionary to keep same order as in the
        # fasta
        from sequana.tools import fast_gc_content

        for chrom, sequence in tqdm(zip(self.fasta.names, self.fasta.sequences)):
            GCs.append(fast_gc_content(sequence))
            chrom_names.append(chrom)
        return chrom_names, GCs

    def plot_GC_per_chromosome(self):
        chrom_names, GCs = self.get_gc_per_chromosome()
        pylab.plot(range(0, len(GCs)), GCs, "o-")
        pylab.xlabel("chromosome")
        pylab.ylabel("GC")

    @functools.lru_cache(maxsize=None)  # None means no limit on cache size
    def get_all_kmer_counts(self, k=6, reverse=False):
        """Get all kmers from the entire genome"""
        counts = defaultdict(int)
        for sequence in tqdm(self.fasta.sequences):
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
                sequence = sequences[index_seq][row.start : row.stop]
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

        # counts = counts[counts > 20]

        odds = []
        for kmer in counts.index:
            count = counts[kmer]

            A = count / Ngenes
            B = kmers[kmer] / genome_size
            odds.append([kmer, count, kmers[kmer], A, B, A / B])

        odds = pd.DataFrame(odds)
        odds.columns = ["kmer", "count", "count_genome", "ratio1", "ratio2", "odds"]
        return odds

    def find_kmers(self, sequence, pattern=r"GCC[AG]CC"):
        """
        >>> k.find_kmers("GCCACC")
        True
        >>> k.find_kmers("AAAAAA")
        False
        """
        kmers = self.builddata()

        import re

        kozak_pattern = re.compile(pattern)
        match = kozak_pattern.search(sequence)
        if match:
            return True
        else:
            return False

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

    def plot_logo(self, df):

        try:
            Nl = len(df["kozak_left"].iloc[0])
            Nr = len(df["kozak_right"].iloc[0])
        except KeyError:
            N = self.left_kozak
            return {"status": "Warning", "msg": "No data found"}

        pos = [defaultdict(int) for x in range(Nl + Nr)]

        for seql, seqr in zip(df.kozak_left.values, df.kozak_right.values):

            for i in range(Nl):
                pos[i][seql[i]] += 1

            for i in range(Nr):
                pos[Nl + i][seqr[i]] += 1

        logo_data = pd.DataFrame(pos).fillna(0)
        if "N" in logo_data.columns:
            del logo_data["N"]

        dd = logo_data.divide(logo_data.sum(axis=1), axis=0)
        self._plot_logo(dd)
        return dd

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

    def _plot_logo(self, data):

        import logomaker

        logo = logomaker.Logo(data)
        for x in [0.25, 0.5, 0.75]:
            pylab.axhline(x, color="grey", zorder=-1, ls="--")
        pylab.axvline(self._left_kozak - 0.5, color="k", lw=2)
        # , color_scheme={'A': 'blue',
        #     'C': 'yellow',
        #      'G': 'green',
        #      'T': 'red'}
        # )

    def plot_kozak_chi2(self, df, GC=0.5):

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

        df = observed
        fig, ax = pylab.subplots(figsize=(8, 6))
        # Plotting values
        ax.matshow(df[["A", "C", "G", "T"]].T, cmap="viridis", vmin=0, vmax=50)

        # Adding text annotations
        df = df[["A", "C", "G", "T"]].T
        for i in range(df.shape[0]):
            for j in range(df.shape[1]):
                ax.text(j, i, f"{df.values[i, j]:.2f}", ha="center", va="center", color="red", fontsize=6)
        pylab.yticks([0, 1, 2, 3], ["A", "C", "G", "T"])
        pylab.xticks(range(len(observed)), ["**" if x < 0.05 else "" for x in observed["chi2"]], fontsize=6)
        pylab.axvline(self.left_kozak - 0.5, lw=2, color="k")

        return observed
