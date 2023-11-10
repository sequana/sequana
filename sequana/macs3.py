#
#  This file is part of Sequana software
#
#  Copyright (c) 2016-2023 - Sequana Development Team
#
#  Distributed under the terms of the 3-clause BSD license.
#  The full license is in the LICENSE file, distributed with this software.
#
#  website: https://github.com/sequana/sequana
#  documentation: http://sequana.readthedocs.io
#
##############################################################################

import colorlog

from sequana.lazy import pandas as pd
from sequana.lazy import pylab
from sequana.utils.pandas import PandasReader
from sequana.viz.venn import plot_venn

logger = colorlog.getLogger(__name__)


__all__ = ["MACS3Reader", "PeakConsensus"]


class MACS3Reader:
    """This class reads output of macs3 tool

    ::

        from sequana import MACS3Reader
        mr = MACS3Reader(data)
        mr.df
        mr.plot_volcano()


    If input file ends in .xls, we assume this is a NAME_peaks.xls file. It is
    therefore a tabulated file where columns are:

    - chromosome name
    - start position of peak
    - end position of peak
    - length of peak region
    - absolute peak summit position
    - pileup height at peak summit
    - -log10(pvalue) for the peak summit (e.g. pvalue =1e-10, then this value should be 10)
    - fold enrichment for this peak summit against random Poisson distribution with local lambda,
    - -log10(qvalue) at peak summit

    Note that coordinates in XLS is 1-based which is different from BED format
    that follows. When --broad is enabled for broad peak calling, the pileup, p-value,
    q-value, and fold change in the XLS file will be the mean value across
    the entire peak region, since peak summit won't be called in broad peak calling mode.

    NAME_peaks.narrowPeak is BED6+4 format file which contains the peak locations
    together with peak summit, p-value, and q-value.

    The fifth column is the **score** calculate as int(-10*log10_pvalue) or
    int(-10*log10_pvalue) where log10-pvalue/log10-qvalue is the 8th or 9th column
    depending on whether -p (pvalue) or -q (qvalue) is used as score cutoff. The 7th
    columns stores the fold-change at peak summit. The 8th is -log10pvalue at peak summit
    and the 9th is -log10qvalue at peak summit. The 10th is the relative summit position
    to peak start equivalent to absolute peak summit position in the xls file.

    NAME_peaks.broadPeak is in BED6+3 format which is similar to narrowPeak file, except
    for missing the 10th column for annotating peak summits. This file and the
    gappedPeak file will only be available when --broad is enabled. Since in the
    broad peak calling mode, the peak summit won't be called, the values in the 5th,
    and 7-9th columns are the mean value across all positions in the peak region.
    Refer to narrowPeak if you want to fix the value issue in the 5th column

    .. reference:: https://github.com/macs3-project/MACS/blob/master/docs/callpeak.md
    """

    def __init__(self, filename):
        # .xls is not stuctured similarly to .narrowPeak
        if filename.endswith(".xls"):
            #
            self.df = PandasReader(filename, comment="#", sep="\t").df
            # narrow and broad files have this header in common and can be checked.
            for x in [
                "chr",
                "start",
                "end",
                "abs_summit",
                "-log10(pvalue)",
                "fold_enrichment",
                "-log10(qvalue)",
                "name",
            ]:
                # we could have also pileup depending on the user's options
                assert x in self.df.columns
            self.df["stop"] = self.df["end"]
        elif filename.endswith("narrowPeak"):
            columns = [
                "chr",
                "start",
                "stop",
                "name",
                "score",
                "NA",
                "fold_enrichment",
                "-log10(pvalue)",
                "-log10(qvalue)",
                "abs_summit_from_start",
            ]
            self.df = PandasReader(filename, header=None, sep="\t", columns=columns).df
        elif filename.endswith("broadPeak"):
            columns = [
                "chr",
                "start",
                "stop",
                "name",
                "score",
                "NA",
                "fold_enrichment",
                "-log10(pvalue)",
                "-log10(qvalue)",
            ]
            self.df = PandasReader(filename, header=None, sep="\t", columns=columns).df
        self.df["length"] = self.df["stop"] - self.df["start"]
        #

    def __len__(self):
        return len(self.df)

    def plot_volcano(self, plotly=False, marker_color="b", title=""):
        """Plots -log10 p-value versus fold change of all peaks"""

        if plotly:
            from plotly import express as px

            df = self.df.copy()
            df["log_adj_pvalue"] = self.df["-log10(pvalue)"]
            df["log2FoldChange"] = pylab.log2(self.df["fold_enrichment"])
            df["info"] = df["chr"] + ":" + df["start"].astype(str) + "-" + df["stop"].astype(str)

            fig = px.scatter(
                df,
                x="log2FoldChange",
                y="log_adj_pvalue",
                hover_name="info",
                log_y=False,
                color="length",
                height=600,
                title=title,
                labels={"log2_fold_enrichment": "log10 p-value"},
            )
            return fig

        pylab.scatter(
            pylab.log2(self.df["fold_enrichment"]),
            self.df["-log10(pvalue)"],
            marker="o",
            alpha=0.5,
            color=marker_color,
            lw=0,
            s=self.df["length"] / self.df["length"].max() * 400,
        )
        pylab.xlabel("Fold enrichment")
        pylab.ylabel("log10 pvalue")
        pylab.title(title)


class PeakConsensus:
    def __init__(self, f1, f2):
        self.df1 = MACS3Reader(f1).df
        self.df1["category"] = "first"
        self.df2 = MACS3Reader(f2).df
        self.df2["category"] = "second"
        self.df_merged = self.merge()

    def merge(self, overlap=0.2):
        df = pd.concat([self.df1, self.df2]).sort_values(["chr", "start"])
        # if overlap at least one base, we merge the peaks and label them with
        # common information, otherwise we report the original peak

        merged = []
        prev = None
        N1 = 0
        N2 = 0
        N12 = 0
        skip_next = True
        for k, current in df.iterrows():
            if skip_next:
                prev = current
                skip_next = False
                continue

            # if current overlaps the prev start or end, there is overlap
            # or if current included in prev there current and prev overlaps
            if current["start"] <= prev["start"] and current["stop"] >= prev["start"]:
                overlap = True
                N12 += 1
            elif current["start"] <= prev["stop"] and current["stop"] >= prev["stop"]:
                overlap = True
                N12 += 1
            elif current["start"] >= prev["start"] and current["stop"] <= prev["stop"]:
                overlap = True
                N12 += 1
            else:
                overlap = False
                if prev["name"].startswith("1_vs_6_7"):
                    N1 += 1
                elif prev["name"].startswith("2_vs_6_7"):
                    N2 += 1

            if overlap:
                m = min(current["start"], prev["start"])
                M = max(current["stop"], prev["stop"])
                data = current.copy()
                data["start"] = m
                data["stop"] = M
                data["category"] = "both"
                merged.append(data)
                skip_next = True
            else:
                m = min(current["start"], prev["start"])
                M = max(current["stop"], prev["stop"])
                merged.append(prev)
                skip_next = False

            prev = current
        df = pd.DataFrame(merged)
        df = df.reset_index(drop=True)
        return df

    def plot_venn(self, title="", labels=[]):
        plot_venn(
            (
                set(self.df_merged.query("category in ['first', 'both']").index),
                set(self.df_merged.query("category in ['second', 'both']").index),
            ),
            labels=(("first", "second")),
        )

    def to_saf(self, filename):
        """For now, all strand are categorised as strand + .
        strand not used in the pipeline for now."""
        df = self.df_merged.reset_index()
        df["strand"] = "+"
        df["GeneID"] = df["index"]
        df[["GeneID", "chr", "start", "stop", "strand", "category"]].to_csv(filename, sep="\t", index=None)

    def to_bed(self, filename):
        # output to be used by homer
        # annotatePeaks.pl file.bed file.fa -gid -gff file.gff
        # GeneID seems to have to start with Interval_ ??

        df = self.df_merged.reset_index()
        df["strand"] = "+"
        df["GeneID"] = ["Interval_{}".format(x) for x in df["index"]]
        mapper = {"both": 3, "first": 1, "second": 2}
        df["score"] = [mapper[x] for x in df["category"]]
        df[["chr", "start", "stop", "GeneID", "score", "strand"]].to_csv(filename, sep="\t", index=None, header=False)
