# -*- coding: utf-8 -*-
#
#  This file is part of Sequana software
#
#  Copyright (c) 2016 - Sequana Development Team
#
#  File author(s):
#      Thomas Cokelaer <thomas.cokelaer@pasteur.fr>
#      Etienne Kornobis <etienne.kornobis@pasteur.fr>
#
# This is a copy of the biokit.viz.volcano module (from myself) since
# biokit will not be maintained in the future (merging its content
# into sequana.
#
#
#  Distributed under the terms of the 3-clause BSD license.
#  The full license is in the LICENSE file, distributed with this software.
#
#  website: https://github.com/sequana/sequana
#  documentation: http://sequana.readthedocs.io
#
##############################################################################
"""Volcano plot"""

from sequana.lazy import pandas as pd
from sequana.lazy import pylab as pylab
from sequana.lazy import numpy as np

from adjustText import adjust_text

__all__ = ["Volcano"]


class Volcano(object):
    """Volcano plot

    In essence, just a scatter plot with annotations.

    .. plot::
        :width: 80%
        :include-source:

        import numpy as np
        fc = np.random.randn(1000)
        pvalue = np.random.randn(1000)
        import pandas as pd
        data = pd.DataFrame([fc, pvalue])
        data = data.T
        data.columns = ['log2FoldChange', 'padj']

        from sequana.viz import Volcano
        v = Volcano(data)
        v.plot()


    """

    def __init__(
        self,
        data=None,
        log2fc_col="log2FoldChange",
        pvalues_col="padj",
        annot_col="",
        color="auto",
        pvalue_threshold=-np.log10(0.05),
        log2fc_threshold=1,
    ):
        """.. rubric:: constructor

        :param DataFrame data: Pandas DataFrame with rnadiff results.
        :param log2fc_col: Name of the column with log2 Fold changes.
        :param pvalues_col:  Name of the column with adjusted pvalues.
        :param annot_col: Name of the column with genes names for plot annotation.
        :param color: for color choice
        :param pvalue_threshold: Adjusted pvalue threshold to use for coloring/annotation.
        :param log2fc_threshold: Log2 Fold Change threshold to use for coloring/annotation.

        """

        self.color = color
        self.pvalue_threshold = pvalue_threshold
        self.log2fc_threshold = log2fc_threshold

        df = data.copy()
        df.rename(
            columns={
                log2fc_col: "log2fc",
                annot_col: "annot",
            },
            inplace=True,
        )
        df["log10pval"] = -np.log10(df[pvalues_col])
        self.df = df
        # self.df = df.loc[df.log10pval.notna()]

        self._get_colors()

    def _get_colors(self):
        """Add colors according to pvalue and fold change thresholds"""

        def coloring(row):
            if (row.log2fc <= -self.log2fc_threshold) and (
                row.log10pval >= self.pvalue_threshold
            ):
                return "royalblue"
            if (row.log2fc >= self.log2fc_threshold) and (
                row.log10pval >= self.pvalue_threshold
            ):
                return "firebrick"
            if (abs(row.log2fc) < self.log2fc_threshold) and (
                row.log10pval >= self.pvalue_threshold
            ):
                return "black"
            else:
                return "lightgrey"

        if self.color == "auto":
            self.df["color"] = self.df.apply(coloring, axis=1)
        else:  # pragma: no cover
            self.df["color"] = self.color

    def annotate(self, **kwargs):
        texts = []
        for row in self.df.itertuples(index=False):
            if not row.annot:
                continue

            if (row.log10pval >= self.pvalue_threshold) and (
                abs(row.log2fc) >= self.log2fc_threshold
            ):
                texts.append(pylab.text(row.log2fc, row.log10pval, row.annot))

        adjust_text(texts, **kwargs)

    def plot(
        self,
        size=10,
        alpha=0.7,
        marker="o",
        fontsize=16,
        xlabel="fold change",
        logy=False,
        threshold_lines={"color": "black", "ls": "--", "width": 0.5},
        ylabel="p-value",
        add_broken_axes=False,
        broken_axes={"ylims": ((0, 10), (50, 100))},
    ):
        """

        :param size: size of the markers
        :param alpha: transparency of the marker
        :param fontsize:
        :param xlabel:
        :param ylabel:
        :param center: If centering the x axis
        """
        pylab.clf()

        if add_broken_axes:  # pragma: no cover
            from brokenaxes import brokenaxes

            _ylims = broken_axes.get("ylims", None)
            _xlims = broken_axes.get("xlims", None)
            bax = brokenaxes(ylims=_ylims, xlims=_xlims)
        else:
            bax = pylab

        bax.scatter(
            self.df.log2fc,
            self.df.log10pval,
            s=size,
            alpha=alpha,
            c=self.df.color,
            marker=marker,
            edgecolors="None",
        )

        # bax.grid()
        # pylab.ylim([0, pylab.ylim()[1]])
        # M = max(abs(self.log2fc)) * 1.1
        # pylab.xlim([-M, M])
        try:
            bax.set_xlabel(xlabel, fontsize=fontsize)
            bax.set_ylabel(ylabel, fontsize=fontsize)
        except:
            bax.xlabel(xlabel, fontsize=fontsize)
            bax.ylabel(ylabel, fontsize=fontsize)

        bax.axhline(
            self.pvalue_threshold,
            color=threshold_lines["color"],
            linestyle=threshold_lines["ls"],
            linewidth=threshold_lines["width"],
        )
        bax.axvline(
            self.log2fc_threshold,
            color=threshold_lines["color"],
            linestyle=threshold_lines["ls"],
            linewidth=threshold_lines["width"],
        )
        bax.axvline(
            -1 * self.log2fc_threshold,
            color=threshold_lines["color"],
            linestyle=threshold_lines["ls"],
            linewidth=threshold_lines["width"],
        )

        if logy is True:
            ax = pylab.gca()
            ax.set(yscale="log")
