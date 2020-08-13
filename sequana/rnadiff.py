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

import glob


__all__ = ["RNADiffResults", "RNADiffResultsOld"]


class RNADiffResults():
    """An object representation of results coming from a RNADiff analysis."""

    def __init__(self, filename, alpha=0.05, log2_fc=0, pattern="*complete*xls",
        sep="\t"):
        """.. rubric:: constructor

        :param rnadiff_results: can be a folder for a simple comparison, or an
            output file containing the rseults of a rnadiff analysis
        :param alpha:
        :param out_dir:
        :param log2_fc: the log2 fold change to apply

        """
        if os.path.isdir(filename):
            filenames = glob.glob(filename + "/tables/" + pattern)
            if len(filenames) == 1:
                self.filename = filenames[0]
            else:
                raise IOError("Found several file with the {} pattern. Please be"
                    "more restrictive using the pattern argument")
        elif os.path.exists(filename):
            # must be a 'complete' file from rnadiff analysis
            self.filename = filename
        else:
            raise TypeError("{} does not exists".format(filename))

        # Finally, read the data itself
        self.df = pd.read_csv(self.filename, sep, index_col=0)

        # Just an alias to a subset of the dataframe
        normed = [x for x in self.df.columns if x.startswith('norm')]
        self.normcounts = self.df[normed]

        # some parameters/attributes
        self._alpha = alpha
        self._log2_fc = log2_fc

        # What are the sample names and condition names
        self.sample_names = [x.replace('norm.', '') for x in normed]
        self.condition_names = set([x[0:-1] for x in self.sample_names])
        self.set_colors()

        self._set_gene_lists_one_condition()

    def set_colors(self, colors=None):
        self.colors = {}
        if colors is None:
            colors = ['orange', 'b', 'r', "y", "k"]
        for i,name in enumerate(self.condition_names):
            try:
                self.colors[name] = colors[i]
            except:
                self.colors[name] = colors[len(colors)-1]

    def _set_log2_fc(self, value):
        self._log2_fc = value
        self._set_gene_lists_one_condition()
    def _get_log2_fc(self):
        return self._log2_fc
    log2_fc = property(_get_log2_fc, _set_log2_fc)

    def _set_alpha(self, value):
        self._alpha = value
    def _get_alpha(self):
        return self._alpha
    alpha = property(_get_alpha, _set_alpha)

    def _set_gene_lists_one_condition(self):

        gene_sets =  {}
        key = "vs".join(sorted(self.condition_names))
        gene_sets = {}

        condition_up = np.logical_and(self.df.log2FoldChange > self.log2_fc,
            self.df.padj <self.alpha)
        gene_sets['up'] = set(self.df[condition_up].index)

        condition_down = np.logical_and(self.df.log2FoldChange < -self.log2_fc,
            self.df.padj <self.alpha)
        gene_sets['down'] = set(self.df[condition_down].index)

        condition_all = np.logical_or(
            np.logical_and(self.df.padj < self.alpha, self.df.log2FoldChange > self.log2_fc),
            np.logical_and(self.df.padj < self.alpha, self.df.log2FoldChange < -self.log2_fc)
        )
        gene_sets['all'] = set(self.df[condition_all].index)

        self.gene_lists = gene_sets

    def summary(self):
        """ Get a summary DataFrame from a RNADiff analysis.
        """
        summary = pd.DataFrame(
            {
                (x,len(self.gene_lists[x])) for x in self.gene_lists.keys()
            }
        )

        df = summary.set_index(0)
        df.columns = ["_vs_".join(self.condition_names)]
        return df

    def plot_count_per_sample(self, fontsize=12, sample_list=None):
        """"Number of mapped reads per sample. Each color for each replicate

        .. plot::
            :include-source:

            from sequana.rnadiff import RNADiffResults
            from sequana import sequana_data

            r = RNADiffResults(sequana_data("rnadiff/rnadiff_onecond_1"))
            r.plot_count_per_sample()
        """
        sample_names = self.sample_names
        N = len(sample_names)
        dd = self.df[sample_names].sum()
        pylab.clf()

        colors = []
        for sample in self.sample_names:
            colors.append(self.colors[self.get_cond_from_sample(sample)])

        pylab.bar(range(N), (dd/1000000).values, 
            color=colors, alpha=1, 
            zorder=10, lw=1, ec="k", width=0.9)
        pylab.xlabel("Samples", fontsize=fontsize)
        pylab.ylabel("Total read count (millions)", fontsize=fontsize)
        pylab.grid(True, zorder=0)
        pylab.title("Total read count per sample", fontsize=fontsize)
        pylab.xticks(range(N), self.sample_names)

    def plot_percentage_null_read_counts(self):
        """

        Bars represent the percentage of null counts in each samples. 
        The dashed horizontal line represents the percentage of 
        feature counts being equal to zero across all samples.

        .. plot::
            :include-source:
    
            from sequana.rnadiff import RNADiffResults
            from sequana import sequana_data

            r = RNADiffResults(sequana_data("rnadiff/rnadiff_onecond_1"))
            r.plot_percentage_null_read_counts()


        """
        N = len(self.sample_names)

        data = (self.df[self.sample_names]==0).sum() 
        data = data / len(self.df) * 100

        all_null = (self.df[self.sample_names].sum(axis=1) == 0).sum()

        colors = []
        for sample in self.sample_names:
            colors.append(self.colors[self.get_cond_from_sample(sample)])

        pylab.clf()
        pylab.bar(range(N), data, 
            color=colors, alpha=1, 
            zorder=10, lw=1, ec="k", width=0.9)
        pylab.axhline(all_null / len(self.df) * 100, lw=2, ls="--", color="k",
            zorder=20)
        pylab.xticks(range(N), self.sample_names)
        pylab.xlabel("Sample")
        pylab.ylabel("Proportion of null counts (%)")
        pylab.grid(True, zorder=0)


    def plot_volcano(self, padj=0.05):
        """


        .. plot::
            :include-source:

            from sequana.rnadiff import RNADiffResults
            from sequana import sequana_data

            r = RNADiffResults(sequana_data("rnadiff/rnadiff_onecond_1"))
            r.plot_volcano()

        """
        d1 = self.df.query("padj>@padj")
        d2 = self.df.query("padj<=@padj")

        fig = pylab.figure()
        pylab.plot(d1.log2FoldChange, -np.log10(d1.padj), marker="o",
            alpha=0.5, color="r", lw=0)
        pylab.plot(d2.log2FoldChange, -np.log10(d2.padj), marker="o",
            alpha=0.5, color="k", lw=0)

        pylab.grid(True)
        pylab.xlabel("fold change")
        pylab.ylabel("log10 adjusted p-value")
        m1 = abs(min(self.df.log2FoldChange))
        m2 = max(self.df.log2FoldChange)
        limit = max(m1,m2)
        pylab.xlim([-limit, limit])
        y1,y2 = pylab.ylim()
        pylab.ylim([0,y2])

        pylab.axhline(-np.log10(0.05), lw=2, ls="--", color="r", label="pvalue threshold (0.05)")

    def plot_pca(self, n_components=2, colors=None):
        """

        .. plot::
            :include-source:

            from sequana.rnadiff import RNADiffResults
            from sequana import sequana_data

            path = sequana_data("rnadiff/rnadiff_onecond_1")
            r = RNADiffResults(path)

            colors = {
                'surexp1': 'r',
                'surexp2':'r',
                'surexp3':'r',
                'surexp1': 'b',
                'surexp2':'b', 
                'surexp3':'b'}
            r.plot_pca(colors=colors)

        """
        from sequana.viz import PCA
        p = PCA(self.df[self.sample_names])
        if colors is None:
            colors = {}
            for sample in self.sample_names:
                colors[sample] = self.colors[self.get_cond_from_sample(sample)]
        p.plot(n_components=n_components, colors=colors)

    def plot_mds(self, n_components=2, colors=None, clf=True):
        from sequana.viz.mds import MDS
        p = MDS(self.df[self.sample_names])
        if colors is None:
            colors = {}
            for sample in self.sample_names:
                colors[sample] = self.colors[self.get_cond_from_sample(sample)]
        p.plot(n_components=n_components, colors=colors, clf=clf)

    def plot_isomap(self, n_components=2, colors=None):
        from sequana.viz.isomap import Isomap
        p = Isomap(self.df[self.sample_names])
        if colors is None:
            colors = {}
            for sample in self.sample_names:
                colors[sample] = self.colors[self.get_cond_from_sample(sample)]
        p.plot(n_components=n_components, colors=colors)


    def get_cond_from_sample(self, sample_name):
        try:
            candidates = [x for x in self.condition_names if sample_name.startswith(x)]
            if len(candidates) == 1:
                return candidates[0]
            else:
                raise ValueError("ambiguous sample name found in several conditions")
        except:
            logger.warning("{} not found".format(sample_name))
            return None




class RNADiffResultsOld(object):  #pragma: no cover
    """ An object representation of results coming from a RNADiff analysis.
    """
    TABLE_DIR_NAME = "tables"
    # File name for the complete table is of form:
    # B1234-v1 (Number of project, number of version)
    SEP = "\t"

    def __init__(self, rnadiff_folder, alpha=0.05, out_dir="gsea", fc=0,
        pattern="*complete.xls"):
        """.. rubric:: constructor

        :param rnadiff_folder:
        :param alpha:
        :param out_dir:
        :param fc:
        :param pattern: to replace the default pattern 

        """

        self.analysis_folder = Path(rnadiff_folder)
        self.name = self.analysis_folder.stem
        self._table_folder = self.analysis_folder / self.TABLE_DIR_NAME

        self.df = self._get_table(pattern)

        # Just an alias to a subset of the dataframe
        normed = [x for x in self.df.columns if x.startswith('norm')]
        self.normcounts = self.df[normed]

        # get the different conditions
        self.comparisons = self.get_comparisons()

        # some parameters
        self.out_dir = out_dir
        self.alpha = alpha

        # What are the sample names and condition names
        self.sample_names = [x.replace('norm.', '') for x in normed]
        self.condition_names = set([x[0:-1] for x in self.sample_names])
        self.set_colors()

        # 
        self.dr_gene_lists = self.get_gene_lists(alpha=alpha)
        if len(self.dr_gene_lists) == 0:
            self.dr_gene_lists = self.get_gene_lists_one_condition(fc=fc, alpha=alpha)

    def set_colors(self, colors=None):
        self.colors = {}
        if colors is None:
            colors = ['orange', 'b', 'r', "y", "k"]
        for i,name in enumerate(self.condition_names):
            try:
                self.colors[name] = colors[i]
            except:
                self.colors[name] = colors[len(colors)-1]


    def _get_table(self, pattern):
        """ Extract the complete (with all comparisons) table 
        from a RNADiff analysis or the normCounts table 
        depending on the pattern specified.
        """
        if pattern:
            table_files = [f for f in self._table_folder.glob(pattern)]
            if len(table_files) == 0:
                raise ValueError("Found no file for your pattern {}".format(pattern))
            elif len(table_files) != 1:
                print(table_files)
                raise ValueError("Found more than 1 file with the pattern {}".format(pattern))
            return pd.read_csv(table_files[0], self.SEP, index_col=0)
        else:
            table_files = [f for f in self._table_folder.glob("*.xls")]
            table = [f for f in table_files if re.match(pattern, str(f))]

            if len(table) == 1 and table[0].is_file():
                return pd.read_csv(table[0], self.SEP, index_col=0)
            else:
                raise IOError(
                    f"Cannot find a single proper table with pattern: {pattern} from RNADiff: {table}"
                )

    def get_comparisons(self):
        """ Get a list of comparisons performed in the RNADiff analysis.
        """

        comparisons = set(
            [x.split(".")[0] for x in self.df.filter(regex=(".*vs.*\..*")).columns]
        )
        return comparisons

    def get_gene_lists_one_condition(self, alpha=None, fc=0):
        if alpha is None:
            alpha = 0.05

        gene_sets =  {}
        key = "vs".join(sorted(self.condition_names))
        gene_sets[key] = {}
        
        condition_up = np.logical_and(self.df.log2FoldChange > fc, self.df.padj <alpha)
        gene_sets[key]['up'] = set(self.df[condition_up].index)

        condition_down = np.logical_and(self.df.log2FoldChange < -fc, self.df.padj <alpha)
        gene_sets[key]['down'] = set(self.df[condition_down].index)

        condition_all = self.df.padj <alpha
        gene_sets[key]['all'] = set(self.df[condition_all].index)
        return gene_sets


    def get_gene_lists(self, alpha=0.05, fc=0):
        """ Get differentially regulated gene lists.
        """

        gene_sets = {}

        for compa in self.comparisons:
            gene_sets[compa] = {}
            regex = f"{compa}\..*"

            log2FC_colname = f"{compa}.log2FoldChange"
            padj_colname = f"{compa}.padj"
            sub_df = self.df.filter(regex=regex)

            query_up = f"`{log2FC_colname}` > {fc} and `{padj_colname}` < {alpha}"
            query_down = f"`{log2FC_colname}` < -{fc} and `{padj_colname}` < {alpha}"
            query_all = f"`{padj_colname}` < {alpha}"

            gene_sets[compa]["up"] = set(sub_df.query(query_up).index)
            gene_sets[compa]["down"] = set(sub_df.query(query_down).index)
            gene_sets[compa]["all"] = set(sub_df.query(query_all).index)

        return gene_sets

    def run_enrichr(self, gene_sets, top_term=20):
        """Run enrichr using gseapy.

        gene_sets: A list of databases coming from the database available with
        gseapy (see gseapy.get_library_name(database="Mouse")) top_term: The max
        number of term to plot.
        """

        out_dir_enr = os.path.join(self.out_dir, "enrichr")
        os.makedirs(out_dir_enr, exist_ok=True)
        failed_comparisons = []
        no_dr_genes = []

        for compa in self.dr_gene_lists:
            for direction in self.dr_gene_lists[compa]:

                ensembl_ids = self.dr_gene_lists[compa][direction]
                gene_symbols = [x.upper() for x in self.df.loc[ensembl_ids].gene_name]

                if len(gene_symbols) == 0:
                    no_dr_genes.append(f"{compa}_{direction}")
                    continue

                for gene_set in gene_sets:

                    # print(f"{gene_set} {compa} {direction} {len(gene_symbols)} ")

                    enr = gseapy.enrichr(
                        gene_list=gene_symbols, gene_sets=gene_set, no_plot=True
                    )

                    title = f"{gene_set}_{compa}_{direction}"

                    out_dir_db = os.path.join(out_dir_enr, gene_set)
                    os.makedirs(out_dir_db, exist_ok=True)

                    # gseapy.plot.barplot(enr, cutoff=0.05)
                    try:
                        gseapy.plot.dotplot(
                            enr.res2d,
                            cutoff=self.alpha,
                            top_term=top_term,
                            ofname=os.path.join(out_dir_db, title + ".pdf"),
                            title=title,
                        )

                    except:
                        failed_comparisons.append(f"{gene_set}_{compa}_{direction}")

        print(f"Failed for: {','.join(failed_comparisons)}")
        print(f"No DR genes found for: {','.join(no_dr_genes)}")

    def summary(self):
        """ Get a summary DataFrame from a RNADiff analysis.
        """
        summary = pd.DataFrame(
            {
                k: {x: len(self.dr_gene_lists[k][x]) for x in self.dr_gene_lists[k]}
                for k in self.dr_gene_lists
            }
        )

        return summary

    def venn(self, compa_list, direction="all", prefix=""):
        """
        Plot a venn diagram comparing the list compa_list of dr gene lists.
        compa_list is a list of comparison names from Deseq2 results
        direction specifies either if up/down/all dr genes are considered
        prefix is a string to be added as prefix to the outfile name.

        compa_list can be a list of lists of comparisons to make.
        ie [["WT", "KO1"],["WT", "KO2"]
        """
        from sequana.viz.venn import plot_venn
        # If compa_list is a list of lists of comparison
        if all(isinstance(l, list) for l in compa_list):

            fig, ax = pylab.subplots(6, 1, figsize=(6, 20))
            ax = ax.flat

            for i, c in enumerate(compa_list):

                plot_venn(
                    [self.dr_gene_lists[x][direction] for x in c],
                    [compa_name for compa_name in c],
                    ax=ax[i],
                )
        # If compa is only a list of comparisons
        else:
            plot_venn(
                [self.dr_gene_lists[x][direction] for x in compa_list],
                [compa_name for compa_name in compa_list],
            )
        out_dir = os.path.join(self.out_dir, "vennDiagrams")
        os.makedirs(out_dir, exist_ok=True)
        outfile = os.path.join(out_dir, f"{prefix}vennDiagrams_{direction}.pdf")

        pylab.savefig(outfile, bbox_inches="tight")

    def compare(self, rnadiff_res_obj, make_plot=False, plot_to_file=""):
        """Compare two RNADiffResults objects.  For now this will plot venn
        diagrams for the comparisons considering UP and DOWN regulated
        genes separately

        If make_plot, will plot venn diagrams

        If simple, returns the overlap (regardless of up or down regulated)
        """

        res = {}

        for compa in self.dr_gene_lists:

            # Up
            up1 = self.dr_gene_lists[compa]["up"]
            up2 = rnadiff_res_obj.dr_gene_lists[compa]["up"]

            # Down
            down1 = self.dr_gene_lists[compa]["down"]
            down2 = rnadiff_res_obj.dr_gene_lists[compa]["down"]

            res[compa] = {
                self.name + "_specific_up": self.df.loc[up1 - up2],
                rnadiff_res_obj.name
                + "_specific_up": rnadiff_res_obj.df.loc[up2 - up1],
                self.name + "_specific_down": self.df.loc[down1 - down2],
                rnadiff_res_obj.name
                + "_specific_down": rnadiff_res_obj.df.loc[down2 - down1],
                self.name
                + "_"
                + rnadiff_res_obj.name
                + "_common_up": self.df.loc[up1.intersection(up2)],
                self.name
                + "_"
                + rnadiff_res_obj.name
                + "_common_down": self.df.loc[down1.intersection(down2)],
            }

        res_summary = {}
        for compa in res:
            res_summary[compa] = {}
            for direction in res[compa]:
                res_summary[compa][direction] = res[compa][direction].shape[0]

        print(pd.DataFrame(res_summary).to_string())

    def plot_count_per_sample(self, fontsize=12, sample_list=None):
        """"Number of mapped reads per sample. Each color for each replicate

        .. plot::
            :include-source:
    
            from sequana.rnadiff import RNADiffResults
            from sequana import sequana_data

            r = RNADiffResults(sequana_data("rnadiff/rnadiff_onecond_1"))
            r.plot_count_per_sample()
        """
        sample_names = [x for x in self.df.columns if x.startswith("norm")] 
        sample_names = [x.replace("norm.", "") for x in sample_names]
        N = len(sample_names)
        dd = self.df[sample_names].sum()
        pylab.clf()
        pylab.bar(range(N), (dd/1000000).values, color=['r']*3+['b']*3, alpha=1)
        pylab.xlabel("Samples", fontsize=fontsize)
        pylab.ylabel("Total read count (millions)", fontsize=fontsize)
        pylab.grid(True)
        pylab.title("Total read count per sample", fontsize=fontsize)

    def plot_percentage_null_read_counts(self):
        """


        Bars represent the percentage of null counts in each samples. 
        The dashed horizontal line represents the percentage of 
        feature counts being equal to zero across all samples.

        .. plot::
            :include-source:
    
            from sequana.rnadiff import RNADiffResults
            from sequana import sequana_data

            r = RNADiffResults(sequana_data("rnadiff/rnadiff_onecond_1"))
            r.plot_percentage_null_read_counts()


        """
        N = len(self.sample_names)

        data = (self.df[self.sample_names]==0).sum() 
        data = data / len(self.df) * 100

        all_null = (self.df[self.sample_names].sum(axis=1) == 0).sum()

        pylab.clf()
        pylab.bar(range(N), data)
        pylab.axhline(all_null / len(self.df) * 100, lw=2, ls="--", color="k")
        pylab.xticks(range(N), self.sample_names)
        pylab.xlabel("Sample")


    def plot_volcano(self):
        """
        .. plot::
            :include-source:
    
            from sequana.rnadiff import RNADiffResults
            from sequana import sequana_data

            r = RNADiffResults(sequana_data("rnadiff/rnadiff_onecond_1"))
            r.plot_volcano()

        """
        d1 = self.df.query("padj>0.05")
        d2 = self.df.query("padj<=0.05")

        fig = pylab.figure()
        pylab.plot(d1.log2FoldChange, -np.log10(d1.padj), marker="o",
            alpha=0.5, color="r", lw=0)
        pylab.plot(d2.log2FoldChange, -np.log10(d2.padj), marker="o",
            alpha=0.5, color="k", lw=0)

        pylab.grid(True)
        pylab.xlabel("fold change")
        pylab.ylabel("log10 adjusted p-value")
        m1 = abs(min(self.df.log2FoldChange))
        m2 = max(self.df.log2FoldChange)
        limit = max(m1,m2)
        pylab.xlim([-limit, limit])
        y1,y2 = pylab.ylim()
        pylab.ylim([0,y2])

        pylab.axhline(-np.log10(0.05), lw=2, ls="--", color="r", label="pvalue threshold (0.05)")

    def pca(self, n_components=2, colors=None):
        """

        .. plot::
            :include-source:

            from sequana.rnadiff import RNADiffResults
            from sequana import sequana_data

            path = sequana_data("rnadiff/rnadiff_onecond_1")
            r = RNADiffResults(path)

            colors = {
                'surexp1': 'r',
                'surexp2':'r',
                'surexp3':'r',
                'surexp1': 'b',
                'surexp2':'b', 
                'surexp3':'b'}
            r.pca(colors=colors)

        """
        from sequana.viz import PCA
        p = PCA(self.df[self.sample_names])
        if colors is None:
            colors = {}
            for sample in self.sample_names:
                colors[sample] = self.colors[self.get_cond_from_sample(sample)]
        p.plot(n_components=n_components, colors=colors)


    def get_cond_from_sample(self, sample_name):
        try:
            candidates = [x for x in self.condition_names if sample_name.startswith(x)]
            if len(candidates) == 1:
                return candidates[0]
            else:
                raise ValueError("ambiguous sample name found in several conditions")
        except:
            logger.warning("{} not found".format(sample_name))
            return None

