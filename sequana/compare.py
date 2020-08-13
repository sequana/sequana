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

    def __init__(self, r1, r2, r3=None):

        if isinstance(r1, RNADiffResults):
            self.r1 = r1
        elif os.path.exists(r1):
            self.r1 = RNADiffResults(r1)
        else:
            raise NotImplementedError

        if isinstance(r2, RNADiffResults):
            self.r2 = r2
        elif os.path.exists(r2):
            self.r2 = RNADiffResults(r2)
        else:
            raise NotImplementedError

        if r3 is None:
            self.r3 = None
        elif isinstance(r3, RNADiffResults):
            self.r3 = r3
        elif os.path.exists(r3):
            self.r3 = RNADiffResults(r3)
        else:
            raise NotImplementedError


    def summary(self):
        conditions = self._get_conditions()
        d1 = self.r1.gene_lists['down']
        d2 = self.r2.gene_lists['down']
        u1 = self.r1.gene_lists['up']
        u2 = self.r2.gene_lists['up']
        res = {
            "up1": len(u1), "up2": len(u2),
            "down1": len(d1), "down2": len(d2),
            "common_down_r1_r2": len(set(d1).intersection(set(d2))),
            "common_up_r1_r2": len(set(u1).intersection(set(u2)))
        }

        if len(conditions) == 3:
            d3 = self.r3.gene_lists['down']
            u3 = self.r3.gene_lists['up']

            res['up3'] = len(u3)
            res['down3'] = len(d3)

        
        return res

    def _get_conditions(self, cond1=None, cond2=None):

        cond1 = "_vs_".join(sorted(self.r1.condition_names))
        cond2 = "_vs_".join(sorted(self.r2.condition_names))
        try:
            cond3 = "_vs_".join(sorted(self.r3.condition_names))
            return cond1, cond2, cond3
        except:
            return cond1, cond2

    def venn_down_only(self, cond1=None, cond2=None, cond3=None, labels=None, ax=None,
            title="Down expressed genes", mode="all"):
        kargs = {}
        kargs['title'] = title
        kargs['labels'] = labels
        kargs['cond1'] = cond1
        kargs['cond2'] = cond2
        kargs['cond3'] = cond3
        kargs['ax'] = ax
        kargs['data1'] = self.r1.gene_lists['down']
        kargs['data2'] = self.r2.gene_lists['down']
        if self.r3 and mode=="all":
            kargs['data3'] =  self.r3.gene_lists['down']
        self._venn(**kargs)


    def venn_up_only(self, cond1=None,cond2=None, cond3=None, labels=None,ax=None,
         title="Up expressed genes", mode="all"):
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
        kargs = {}
        kargs['title'] = title
        kargs['labels'] = labels
        kargs['cond1'] = cond1
        kargs['cond2'] = cond2
        kargs['cond3'] = cond3
        kargs['ax'] = ax
        kargs['data1'] = self.r1.gene_lists['up']
        kargs['data2'] = self.r2.gene_lists['up']
        if self.r3 and mode=="all":
            kargs['data3'] =  self.r3.gene_lists['up']
        self._venn(**kargs)

    def _venn(self, data1, data2, data3=None, cond1=None, cond2=None, cond3=None, labels=None,
        ax=None, title="expressed genes"):

        from sequana.viz.venn import plot_venn
        if data3 is None:
            # if mode is two_only, this returns 3 conditions
            #cond1, cond2 = self._get_conditions()
            if labels is None:
                labels = [cond1, cond2]
            plot_venn([data1, data2],
                labels=labels, ax=ax, title=title)
        else:
            cond1, cond2, cond3 = self._get_conditions()
            if labels is None:
                labels = [cond1, cond2, cond3]
            plot_venn(
                [data1, data2, data3],
                labels=labels, ax=ax, title=title)

    def venn_all(self, cond1=None, cond2=None, cond3=None, labels=None, ax=None,
        title="all expressed genes", mode="all"):

        kargs = {}
        kargs['title'] = title
        kargs['labels'] = labels
        kargs['cond1'] = cond1
        kargs['cond2'] = cond2
        kargs['cond3'] = cond3
        kargs['ax'] = ax
        kargs['data1'] = self.r1.gene_lists['all']
        kargs['data2'] = self.r2.gene_lists['all']
        if self.r3 and mode=="all":
            kargs['data3'] =  self.r3.gene_lists['all']
        self._venn(**kargs)

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
        #cond1, cond2 = self._get_cond1_cond2()
        if labels is None:
            labels = ['r1', 'r2']

        if mode in ["down"]:
            # Negative values !
            A = self.r1.df.loc[self.r1.gene_lists[mode]].sort_values(by=sortby)
            B = self.r2.df.loc[self.r2.gene_lists[mode]].sort_values(by=sortby)
        else:
            A = self.r1.df.loc[self.r1.gene_lists[mode]].sort_values(
                by=sortby, ascending=False)
            B = self.r2.df.loc[self.r2.gene_lists[mode]].sort_values(
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
                self.venn_down_only(ax=ax, title=None, labels=labels,
                    mode="two_only")
            elif mode=="up":
                self.venn_up_only(ax=ax, title=None, labels=labels,
                    mode="two_only")
            elif mode=="all":
                self.venn_all(ax=ax, title=None, labels=labels,
                    mode="two_only")


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
        #cond1, cond2 = self._get_cond1_cond2()
        cond1, cond2 = "cond1", "cond2"
        if labels is None:
            #labels = [cond1, cond2]
            labels = [cond1, cond2]

        A = self.r1.df.loc[self.r1.gene_lists["all"]]
        B = self.r2.df.loc[self.r2.gene_lists["all"]]

        if cond1 == cond2:
            cond1 += "(1)"
            cond2 += "(2)"

        pylab.clf()
        pylab.plot(A.log2FoldChange, -np.log10(A.padj), marker="o",
            alpha=0.5, color="r", lw=0, label=labels[0], pickradius=4,
            picker=True)
        pylab.plot(B.log2FoldChange, -np.log10(B.padj), marker="x",
            alpha=0.5, color="k", lw=0, label=labels[1], pickradius=4,
            picker=True)

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



