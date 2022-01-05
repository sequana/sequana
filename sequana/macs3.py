#
#  This file is part of Sequana software
#
#  Copyright (c) 2016-2021 - Sequana Development Team
#
#  Distributed under the terms of the 3-clause BSD license.
#  The full license is in the LICENSE file, distributed with this software.
#
#  website: https://github.com/sequana/sequana
#  documentation: http://sequana.readthedocs.io
#
##############################################################################
import itertools

from sequana.lazy import pandas as pd
from sequana.lazy import pylab
from sequana.lazy import numpy as np

import colorlog

logger = colorlog.getLogger(__name__)


class PandasReader:
    def __init__(self, filename, sep="\t", **kwargs):
        self.df = pd.read_csv(filename, sep=sep, **kwargs)
        # If there is a header, let us strip it from spaces
        try:
            self.df.columns = [x.strip() for x in self.df.columns]
        except:
            pass

        # let us strip strings from spaces
        for x in self.df.columns:
            try:
                self.df[x] = self.df[x].str.strip()
            except Exception:
                pass


class MACS3Reader(PandasReader):
    """This class reads output of macs3 tool


    If input file ends in .xls, we assume this is a NAME_peaks.xls file:
    a tabular file which contains information about called peaks. The columns
    are then filled using the header that should be:

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
    depending on whether -p (pvalue) or -q (qvalue) is used as score cutoff.
    The 7th columns stores the fold-change at peak summit
    8th: -log10pvalue at peak summit
    9th: -log10qvalue at peak summit
    10th: relative summit position to peak start equivalent to absolute peak
    summit position in the xls file

    NAME_peaks.broadPeak is in BED6+3 format which is similar to narrowPeak file, except for missing the 10th column for annotating peak summits. This file and the gappedPeak file will only be available when --broad is enabled. Since in the broad peak calling mode, the peak summit won't be called, the values in the 5th, and 7-9th columns are the mean value across all positions in the peak region. Refer to narrowPeak if you want to fix the value issue in the 5th column
    .. reference:: https://github.com/macs3-project/MACS/blob/master/docs/callpeak.md
    """

    def __init__(self, filename):

        # .xls is not stuctured similarly to .narrowPeak
        if filename.endswith(".xls"):
            #
            self.df = pd.read_csv(filename, comment="#", sep="\t")
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
                # we could have also pileup depeding on the user's options
                assert x in self.df.columns
            self.df["stop"] = self.df["end"]
        elif filename.endswith("narrowPeak"):
            self.df = pd.read_csv(filename, header=None, sep="\t")
            self.df.columns = [
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
            self.df["end"] = self.df["stop"]
        elif filename.endswith("broadPeak"):
            self.df = pd.read_csv(filename, header=None, sep="\t")
            self.df.columns = [
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
            self.df["end"] = self.df["stop"]
        self.df["length"] = self.df["stop"] - self.df["start"]
        #

    def __len__(self):
        return len(self.df)

    def plot_volcano(self, plotly=False, marker_color="b"):

        if plotly:
            from plotly import express as px

            df = self.df.copy()
            df["log_adj_pvalue"] = self.df["-log10(pvalue)"]
            df["log2FoldChange"] = pylab.log2(self.df["fold_enrichment"])
            # df['hover_name'] = df['start']
            hover_name = "start"

            try:
                df["info"] = df["chr"] + ":" + df["start"].astype(str) + "-" + df["stop"].astype(str)
            except KeyError:
                df["info"] = df["chr"] + ":" + df["start"].astype(str) + "-" + df["end"].astype(str)

            fig = px.scatter(
                df,
                x="log2FoldChange",
                y="log_adj_pvalue",
                hover_name="info",
                log_y=False,
                color="length",
                height=600,
                labels={"log2_fold_enrichment": "log10 p-value"},
            )
            """fig.update_layout(
                 shapes=[dict(type='line', 
                               xref='x', x0=df.log2FoldChange.min(), x1=df.log2FoldChange.max(), 
                               yref='y', y0=-pylab.log10(padj), y1=-pylab.log10(padj), 
                             line=dict(
                                  color="black",
                                  width=1,
                                  dash="dash"))
                        ])
            """
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
