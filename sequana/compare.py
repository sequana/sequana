# -*- coding: utf-8 -*-
#
#  This file is part of Sequana software
#
#  Copyright (c) 2016-2020 - Sequana Development Team
#
#  File author(s):
#      Thomas Cokelaer <thomas.cokelaer@pasteur.fr>
#      Etienne Kornobis <etienne.kornobis@pasteur.fr>
#
#  Distributed under the terms of the 3-clause BSD license.
#  The full license is in the LICENSE file, distributed with this software.
#
#  website: https://github.com/sequana/sequana
#  documentation: http://sequana.readthedocs.io
#
##############################################################################

from pathlib import Path
import re
import os

from sequana.lazy import pandas as pd
from sequana.lazy import pylab
from sequana.lazy import numpy as np
from sequana import logger
from matplotlib_venn import venn2_unweighted, venn3_unweighted

from sequana.rnadiff import RNADiffResults


__all__ = ["RNADiffCompare"]



class Compare():
    def __init__(self):
        pass

class RNADiffCompare(Compare):
    """ An object representation of results coming from 
    a RNADiff analysis.

    ::

        from sequana import sequana_data
        from sequana.compare import RNADiffCompare

        c = RNADiffCompare(
            sequana_data("rnadiff/rnadiff_onecond_1"),
            sequana_data("rnadiff/rnadiff_onecond_2"))

        # get number of up/down expressed genes and number of
        # common genes
        print(c.summary())


    """

    def __init__(self, r1, r2):

        if isinstance(r1, RNADiffResults):
            self.r1 = r1
        elif os.path.exists(r1):
            self.r1 = RNADiffResults(r1)

        if isinstance(r2, RNADiffResults):
            self.r2 = r2
        elif os.path.exists(r2):
            self.r2 = RNADiffResults(r2)
        

    def summary(self):
        cond1, cond2 = self._get_cond1_cond2()
        d1 = self.r1.dr_gene_lists[cond1]['down']
        d2 = self.r2.dr_gene_lists[cond2]['down']
        u1 = self.r1.dr_gene_lists[cond1]['up']
        u2 = self.r2.dr_gene_lists[cond2]['up']

        res = {
            "up1": len(u1), "up2": len(u2),
            "down1": len(d1), "down2": len(d2),
            "common_down": len(set(d1).intersection(set(d2))),
            "common_up": len(set(u1).intersection(set(u2)))
        }
        return res

    def _get_cond1_cond2(self, cond1=None, cond2=None):
        k1 = self.r1.dr_gene_lists.keys()
        k2 = self.r2.dr_gene_lists.keys()
        if len(k1) == 1:
            cond1 = list(k1)[0]
        elif cond1: # try it 
            self.r1.dr_gene_lists[cond1]
        else: # pragma: no cover
            raise KeyError("You must set the cond1 name amongst: {}".format(k1))

        if len(k2) == 1:
            cond2 = list(k2)[0]
        elif cond2: # try it 
            self.r2.dr_gene_lists[cond2]
        else:
            raise KeyError("You must set the cond2 name amongst: {}".format(k2))
        return cond1, cond2

    def venn_down_only(self, cond1=None, cond2=None, labels=None, ax=None,
            title="Down expressed genes"):

        cond1, cond2 = self._get_cond1_cond2()
        if labels is None:
            labels = [cond1, cond2]
        from sequana.viz.venn import plot_venn
        plot_venn([
                self.r1.dr_gene_lists[cond1]['down'],
                self.r2.dr_gene_lists[cond2]['down']], 
                labels=labels, ax=ax,
                title=title)

    def venn_up_only(self, cond1=None,cond2=None, labels=None,ax=None,
         title="Up expressed genes"):
        """Venn diagram of cond1 from RNADiff result1 vs cond2 in RNADiff
        result 2

        .. plot::
            :include-source:

            from sequana import sequana_data
            from sequana.compare import RNADiffCompare

            c = RNADiffCompare(
                sequana_data("rnadiff/rnadiff_onecond_1"),
                sequana_data("rnadiff/rnadiff_onecond_2"))
            c.venn_up_only()
        """
        cond1, cond2 = self._get_cond1_cond2()
        if labels is None:
            labels = [cond1, cond2]
        from sequana.viz.venn import plot_venn
        plot_venn([
                self.r1.dr_gene_lists[cond1]['up'],
                self.r2.dr_gene_lists[cond2]['up']], 
                labels=labels, ax=ax, title=title)

    def venn_all(self, cond1=None, cond2=None, labels=None, ax=None, 
        title="all expressed genes"):
        cond1, cond2 = self._get_cond1_cond2()
        if labels is None:
            labels = [cond1, cond2]
        from sequana.viz.venn import plot_venn
        plot_venn([
                self.r1.dr_gene_lists[cond1]['all'],
                self.r2.dr_gene_lists[cond2]['all']], 
                labels=labels, ax=ax, title=title)


    def plot_common_major_counts(self, mode, labels=None,
            switch_up_down_cond2=False, add_venn=True, xmax=None, 
            title="", fontsize=12, sortby="log2FoldChange"):
        """

        :param mode: down, up or all


        .. plot::
            :include-source:

            from sequana import sequana_data
            from sequana.compare import RNADiffCompare

            c = RNADiffCompare(
                sequana_data("rnadiff/rnadiff_onecond_1"),
                sequana_data("rnadiff/rnadiff_onecond_2"))
            c.plot_common_major_counts("down")
        """
        cond1, cond2 = self._get_cond1_cond2()
        if labels is None:
            labels = [cond1, cond2]

        if mode in ["down"]:
            # Negative values !
            A = self.r1.df.loc[self.r1.dr_gene_lists[cond1][mode]].sort_values(by=sortby)
            B = self.r2.df.loc[self.r2.dr_gene_lists[cond2][mode]].sort_values(by=sortby)
        else:
            A = self.r1.df.loc[self.r1.dr_gene_lists[cond1][mode]].sort_values(
                by=sortby, ascending=False)
            B = self.r2.df.loc[self.r2.dr_gene_lists[cond2][mode]].sort_values(
                by=sortby, ascending=False)
        # sometimes, up and down may be inverted as compared to the other
        # conditions
        N = []
        for i in range(1,max(len(A), len(B))):
            a = A.iloc[0:i].index
            b = B.iloc[0:i].index
            n = len(set(b).intersection(set(a)))
            N.append(n / i*100)

        max_common = len(set(A.index).intersection(set(B.index)))
        pylab.clf()
        if len(A) > len(B):
            pylab.axhline(max_common/len(A)*100, color="r", ls='--', label="min set intersection")
            pylab.axvline(len(B), ls="--", color="k", label="rank of minor set")
        else:
            pylab.axhline(max_common/len(B)*100, color='r', ls='--', label="min set intersect")
            pylab.axvline(len(A), ls="--", color="k", label="rank of minor set")

        pylab.plot(N)
        pylab.xlabel('rank', fontsize=fontsize)
        pylab.ylabel('% common features', fontsize=fontsize)
        pylab.grid(True)
        pylab.ylim([0,100])
        if xmax:
            pylab.xlim([0, xmax])
        else:
            pylab.xlim([0, max(len(A),len(B))])
        pylab.title(title, fontsize=fontsize)
        ax = pylab.gca()
        ax2 = ax.twinx()
        ax2.plot(A[sortby].values, "orange", label=sortby)
        ax2.set_ylabel(sortby)
        pylab.legend(loc="lower left")
        ax.legend(loc="lower right")

        if add_venn:
            f = pylab.gcf()
            ax = f.add_axes([0.5,0.5,0.35,0.35], facecolor="grey")
            if mode=="down":
                self.venn_down_only(ax=ax, title=None, labels=labels)
            elif mode=="up":
                self.venn_up_only(ax=ax, title=None, labels=labels)
            elif mode=="all":
                self.venn_all(ax=ax, title=None, labels=labels)


    def plot_volcano(self, labels=None):
        """Volcano plot of log2 fold change versus log10 of adjusted p-value

        .. plot::
            :include-source:

            from sequana import sequana_data
            from sequana.compare import RNADiffCompare

            c = RNADiffCompare(
                sequana_data("rnadiff/rnadiff_onecond_1"),
                sequana_data("rnadiff/rnadiff_onecond_2"))
            c.plot_volcano()
        """
        cond1, cond2 = self._get_cond1_cond2()
        if labels is None:
            labels = [cond1, cond2]
        A = self.r1.df.loc[self.r1.dr_gene_lists[cond1]["all"]]
        B = self.r2.df.loc[self.r2.dr_gene_lists[cond2]["all"]]

        if cond1 == cond2:
            cond1 += "(1)"
            cond2 += "(2)"

        pylab.clf()
        pylab.plot(A.log2FoldChange, -np.log10(A.padj), marker="o",
            alpha=0.5, color="r", lw=0, label=labels[0], picker=4)
        pylab.plot(B.log2FoldChange, -np.log10(B.padj), marker="x",
            alpha=0.5, color="k", lw=0, label=labels[1], picker=4)

        genes = list(A.index) + list(B.index)
        pylab.grid(True)
        pylab.xlabel("fold change")
        pylab.ylabel("log10 adjusted p-value")
        pylab.legend(loc="lower right")
        ax = pylab.gca()

        def onpick(event):
            thisline = event.artist
            self.event = event
            label = thisline.get_label()
            if label == cond1:
                gene_name = A.index[event.ind[0]]
                x1 = round(A.loc[gene_name].log2FoldChange,1)
                y1 = round(-np.log10(A.loc[gene_name].padj),1)
                try:
                    x2 = round(B.loc[gene_name].log2FoldChange,1)
                    y2 = round(-np.log10(B.loc[gene_name].padj),1)
                except:
                    x2, y2 = None, None
            else:
                gene_name = B.index[event.ind[0]]
                x1 = round(B.loc[gene_name].log2FoldChange,1)
                y1 = round(-np.log10(B.loc[gene_name].padj),1)
                try:
                    x2 = round(A.loc[gene_name].log2FoldChange,1)
                    y2 = round(-np.log10(A.loc[gene_name].padj),1)
                except:
                    x2, y2 = None, None

            try:
                if x2 is None:
                    ax.title.set_text("{} at pos [{},{}]".format(
                        gene_name,x1,y1))
                else:
                    ax.title.set_text("{} at pos [{},{}] and [{},{}]".format(
                            gene_name,x1,y1,x2,y2))
            except:
                print("exception")
                ax.title.set_text("")
            pylab.draw()
        fig = pylab.gcf()
        fig.canvas.mpl_connect('pick_event', onpick)



