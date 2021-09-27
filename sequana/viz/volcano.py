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

        from sequana.viz import Volcano
        v = Volcano(fc, -np.log10(pvalue**2))
        v.plot()


    """

    def __init__(
        self,
        fold_changes=None,
        pvalues=None,
        color="auto",
        pvalue_threshold=0.05,
        fold_change_threshold=1,
    ):
        """.. rubric:: constructor


        :param list fold_changes: 1D array or list
        :param list pvalues: 1D array or list
           the threshold provided.
        :param pvalue_threshold: adds an horizontal dashed line at
        :param fold_change_threshold: colors in grey the absolute fold
            changes below a given threshold.
        """

        # try to compute the FC now
        # if self.fold_change is None:
        #    self.fold_change = pylab.log2(X1/X0)

        # if pvalue is None:
        #    # assume a normal distribution mean 0 and sigma 1
        #    import scipy.stats
        #    self.pvalue = - pylab.log10(scipy.stats.norm.pdf(abs(self.fold_change), 0,1)),

        self.fold_changes = np.array(fold_changes)
        self.pvalues = np.array(pvalues)
        self.color = color
        self.pvalue_threshold = pvalue_threshold
        self.fold_change_threshold = fold_change_threshold
        assert len(self.fold_changes) == len(self.pvalues)

        self.df = pd.DataFrame(
            {"fold_change": self.fold_changes, "pvalue": self.pvalues}
        )
        self._get_colors()

    def _get_colors(self):
        """Add colors according to pvalue and fold change thresholds"""

        def coloring(row):
            if (row.fold_change <= -self.fold_change_threshold) and (
                row.pvalue >= self.pvalue_threshold
            ):
                return "royalblue"
            if (row.fold_change >= self.fold_change_threshold) and (
                row.pvalue >= self.pvalue_threshold
            ):
                return "firebrick"
            if (abs(row.fold_change) < self.fold_change_threshold) and (
                row.pvalue >= self.pvalue_threshold
            ):
                return "black"
            else:
                return "lightgrey"

        if self.color == "auto":
            self.df["color"] = self.df.apply(coloring, axis=1)
        else:  # pragma: no cover
            self.df["color"] = self.color

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
            self.df.fold_change,
            self.df.pvalue,
            s=size,
            alpha=alpha,
            c=self.df.color,
            marker=marker,
            edgecolors="None",
        )

        bax.grid()
        # pylab.ylim([0, pylab.ylim()[1]])
        # M = max(abs(self.fold_change)) * 1.1
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
            self.fold_change_threshold,
            color=threshold_lines["color"],
            linestyle=threshold_lines["ls"],
            linewidth=threshold_lines["width"],
        )
        bax.axvline(
            -1 * self.fold_change_threshold,
            color=threshold_lines["color"],
            linestyle=threshold_lines["ls"],
            linewidth=threshold_lines["width"],
        )

        if logy is True:
            ax = pylab.gca()
            ax.set(yscale="log")
