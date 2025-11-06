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
import csv
from pathlib import Path
from typing import List, Optional, Union

import colorlog
import pandas as pd

from sequana.lazy import pylab

logger = colorlog.getLogger(__name__)


__all__ = ["TRF"]


class TRF:  # pragma: no cover
    """Tandem Repeat Finder utilities

    The input data is the output of trf tool when using the -d option.
    This is not a CSV file. It contains comments in the middle of the file to
    indicate the name of the contig.

    The output filename has the following filename convention::

        test.fa.2.5.7.80.10.50.2000.dat

    where the numbers indicate the 7 input parameters:

    * Match  = matching weight
    * Mismatch  = mismatching penalty
    * Delta = indel penalty
    * PM = match probability (whole number)
    * PI = indel probability (whole number)
    * Minscore = minimum alignment score to report
    * MaxPeriod = maximum period size to report

    You may use ``-h`` to suppress html output.

    Then, you can use this class to easly identify the pattern you want::

        t = TRF("input.dat")
        query = "length>100 and period_size==3 and entropy>0 and C>20 and A>20 and G>20"
        t.df.query(query)

    """

    def __init__(self, filename: Union[str, Path], verbose: bool = False, frmt: Optional[str] = None):
        filename = Path(filename)
        self.filename = filename

        if frmt is None:
            if filename.suffix == ".csv":
                frmt = "csv"
            elif filename.suffix == ".dat":
                frmt = "trf"
            else:
                raise ValueError(f"Unknown file type for {filename}. Use frmt='trf' or 'csv'.")
        if frmt == "trf":
            # input can be the output of TRF or our trf dataframe
            self.df = self._parse_trf(verbose=verbose)
        else:
            self.df = pd.read_csv(filename)

    def __repr__(self) -> str:
        N = len(self.df.seq1.unique())
        msg = "Number of unique pattern found: {}\n".format(N)
        msg += "Number of entries: {}".format(len(self.df))
        return msg

    def __len__(self):
        return len(self.df)

    def _parse_trf(self, verbose: bool = True, max_seq_length: int = 20) -> pd.DataFrame:
        """scan output of trf and returns a dataframe

        The format of the output file looks like::

            Tandem Repeats Finder Program

            some info

            Sequence: chr1

            Parameters: 2 5 7 80 10 50 2000

            10001 10468 6 77.2 6 95 3 801 33 51 0 15 1.43 TAACCC TAACCCTA...
            1 10 6 77.2 6 95 3 801 33 51 0 15 1.43 TAACCC TAACCCTA...

            Sequence: chr2

            Parameters: 2 5 7 80 10 50 2000

            10001 10468 6 77.2 6 95 3 801 33 51 0 15 1.43 TAACCC TAACCCTA...

        The dataframe stores a row for each sequence and each pattern found. For
        instance, from the example above you will obtain 3 rows, two for the
        first sequence, and one for the second sequence.
        """
        data = []
        sequence_name = None

        with open(self.filename, "r") as fin:

            # skip lines until we reach "Sequence"
            while sequence_name is None:
                line = fin.readline()
                if line.startswith("Sequence:"):
                    sequence_name = line.split()[1].strip()

            # now we read the rest of the file
            count = 0
            # If we concatenate several files, we also want to ignore the header
            for line in fin.readlines():
                if line.startswith("Sequence:"):
                    sequence_name = line.split()[1].strip()
                    count += 1
                    if count % 100000 == 0:
                        logger.info("scanned {} sequences".format(count))
                    # logger.info("scanned {} sequences".format(count))
                else:
                    this_data = line.split()
                    if len(this_data) == 15:
                        this_data[14] = this_data[14][0:max_seq_length]
                        data.append([sequence_name] + this_data)

        if not data:
            raise ValueError(f"No valid TRF data found in {self.filename}")

        columns = [
            "sequence_name",
            "start",
            "end",
            "period_size",
            "CNV",
            "size_consensus",
            "percent_matches",
            "percent_indels",
            "score",
            "A",
            "C",
            "G",
            "T",
            "entropy",
            "seq1",
            "seq2",
        ]
        df = pd.DataFrame(data, columns=columns)

        numeric_cols = [
            "start",
            "end",
            "period_size",
            "CNV",
            "size_consensus",
            "percent_matches",
            "percent_indels",
            "score",
            "A",
            "C",
            "G",
            "T",
            "entropy",
        ]
        df[numeric_cols] = df[numeric_cols].apply(pd.to_numeric, errors="coerce")

        df["length"] = df["end"] - df["start"] + 1
        return df

    def hist_cnvs(self, bins=50, CNVmin=10, motif=None, color="r", log=True):
        """Plot histogram of CNV counts for a given motif list."""
        motif = motif or ["CAG", "AGC", "GCA"]
        subset = self.df.query("CNV > @CNVmin and seq1 in @motif")
        subset["CNV"].hist(bins=bins, log=log, color=color)
        pylab.xlabel("CNV length (bp)")
        pylab.ylabel("Count")
        pylab.title("Histogram of CNV values")

    def hist_length_repetition(self, bins=50, CNVmin=3, motif=["CAG", "AGC", "GCA"], color="r", log=True):
        """

        histogram of the motif found in the list provided by users.
        As an example, this is triplet CAG. Note that we also add the shifted
        version AGC and GCA.

        """
        cnvs = self.df.query("CNV>@CNVmin and seq1 in @motif").CNV
        repet = self.df.query("CNV>@CNVmin and seq1 in @motif").period_size
        data = cnvs * repet
        data.hist(bins=bins, log=log, color=color)
        pylab.xlabel("Repetition length (bp)")
        pylab.ylabel("#")

    def hist_period_size(self, bins=50):
        """Length of the repetitions"""
        self.df.period_size.hist(bins=bins)
        pylab.xlabel("repeat length")

    def hist_entropy(self, bins=50):
        """Histogram of the entropy of all found repeats"""
        self.df.entropy.hist(bins=bins)
        pylab.xlabel("Entropy")
        pylab.ylabel("#")

    def hist_repetitions_per_sequence(self):
        # How many repetitions per sequence
        counts = self.df.groupby("sequence_name").size()
        counts.hist()
        pylab.xlabel("# repetitions per sequence")
        pylab.ylabel("Count")
        pylab.title("Histogram of repeats per sequence")

    def to_bed(self, outfile, cmap="autumn"):
        import matplotlib.colors as colors
        from matplotlib.cm import get_cmap

        cmap = get_cmap(cmap)
        norm = colors.Normalize(vmin=0, vmax=self.df.length.median() * 2)

        with open(outfile, "w") as fout:
            for _, row in self.df.iterrows():

                length = row["length"]
                rgba = cmap(norm(length))
                CNV = row["CNV"]
                R, G, B = int(rgba[0] * 255), int(rgba[1] * 255), int(rgba[2] * 255)

                data = "\t".join(
                    [
                        row["sequence_name"],
                        str(row["start"]),
                        str(row["end"]),
                        f"CNV={CNV}_length={length}",
                        str(row["score"]),
                        "+",
                        str(row["start"]),
                        str(row["end"]),
                        f"{R},{G},{B}",
                    ]
                )
                fout.write(data + "\n")

    def summary(self):
        """Return a descriptive summary of TRF results."""
        summary = {
            "n_entries": len(self.df),
            "n_sequences": self.df["sequence_name"].nunique(),
            "median_period": self.df["period_size"].median(),
            "median_length": self.df["length"].median(),
            "max_CNV": self.df["CNV"].max(),
            "motif_counts": self.df["seq1"].value_counts().head(10).to_dict(),
        }
        return summary

    def top_motifs(self, n=10):
        """Return the most common repeat motifs."""
        return self.df["seq1"].value_counts().head(n)

    def _get_counts_and_length(self):
        counts = self.df.groupby("sequence_name").size()
        lengths = self.df.groupby("sequence_name")["end"].max()
        return counts, lengths

    def plot_repeat_density(self, normalize_by_length=True):
        """Plot repeat density per sequence."""
        counts, lengths = self._get_counts_and_length()
        if normalize_by_length and "end" in self.df:
            counts = counts / (lengths / 1e3)
            ylabel = "Repeats per kb"
        else:
            ylabel = "Repeats per sequence"
        counts.sort_values(ascending=False).plot(kind="bar", figsize=(10, 4))
        pylab.ylabel(ylabel)
        pylab.xlabel("Sequence")
        pylab.title("Repeat density per sequence")
        pylab.tight_layout()
        return counts, lengths

    def plot_period_vs_CNV(self, log=True):
        """Scatter plot: period size vs copy number."""
        ax = self.df.plot.scatter("period_size", "CNV", alpha=0.4)
        if log:
            ax.set_xscale("log")
            ax.set_yscale("log")
        pylab.title("Period size vs CNV")
        pylab.tight_layout()

    def entropy_distribution_by_period(self):
        """Boxplot of entropy per period size."""
        self.df.boxplot(column="entropy", by="period_size", figsize=(8, 4))
        pylab.suptitle("")
        pylab.title("Entropy distribution per period size")
        pylab.ylabel("Entropy")
        pylab.tight_layout()

    def to_bedgraph(self, outfile):
        """Export TRF data as a bedGraph of repeat density."""
        bg = self.df[["sequence_name", "start", "end", "CNV"]].copy()
        bg.columns = ["chrom", "start", "end", "score"]
        bg.to_csv(outfile, sep="\t", index=False, header=False)
        logger.info(f"Saved bedGraph to {outfile}")

    def get_repeats_in_region(self, seq_name, start, end):
        """Return repeats overlapping a genomic region."""
        mask = (self.df["sequence_name"] == seq_name) & (self.df["end"] >= start) & (self.df["start"] <= end)
        return self.df.loc[mask]

    def plot_motif_logo(self, motif, n=50):
        import logomaker

        """Create sequence logo for a specific motif."""
        seqs = self.df.query("seq1 == @motif")["seq2"].head(n)
        if seqs.empty:
            print("No such motif found.")
            return
        df_logo = logomaker.alignment_to_matrix(seqs.tolist())
        logomaker.Logo(df_logo)
        pylab.title(f"Motif logo for {motif}")

    def plot_count_versus_length(self):
        counts, lengths = self._get_counts_and_length()
        import numpy as np

        X = np.array(counts / lengths * 1000)
        Y = np.array(lengths.values)
        mask = (X > 0) & (Y > 0)
        X, Y = X[mask], Y[mask]

        # Compute regression in log space
        logX = np.log10(X)
        logY = np.log10(Y)
        slope, intercept = np.polyfit(logX, logY, 1)

        # Predicted line in original scale
        Y_pred = 10 ** (intercept + slope * logX)

        # Plot
        pylab.clf()
        pylab.scatter(X, Y, alpha=0.6, label="Data")
        pylab.plot(X, Y_pred, color="red")  # , label=f'Regression: y = {10**intercept:.2f} x^{slope:.2f}')
        pylab.xscale("log")
        pylab.yscale("log")
        pylab.xlabel("Repeats per kb")
        pylab.ylabel("Chromosome size")
        pylab.grid(True, which="both", ls="--", alpha=0.5)
        # from scipy import pearsonr
        #  pearsonr(counts, lengths)
