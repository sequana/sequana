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
from jinja2 import Environment, PackageLoader
import subprocess

from sequana.lazy import pandas as pd
from sequana.lazy import pylab
from sequana.lazy import numpy as np
from sequana import logger
from sequana.gff3 import GFF3

import glob

__all__ = ["RNADiffAnalysis", "RNADiffResults"]


class RNADiffAnalysis:
    """A tool to prepare and run a RNA-seq differential analysis with DESeq2

    :param counts_tsv: Path to tsv file out of FeatureCount with all samples together.
    :param groups_tsv: Path to tsv file with the definition of the groups for each sample.
    :param condition: The name of the column from groups_tsv to use as condition. For more
        advanced design, a R function of the type 'condition*inter' (without the '~') could
        be specified (not tested yet). Each name in this function should refer to column
        names in groups_tsv.
    :param comparisons: A list of tupples indicating comparisons to be made e.g A vs B would be [("A", "B")]
    :param batch: None for no batch effect or name of a column in groups_tsv to add a batch effec.
    :param fit_type: Default "parametric".
    :param beta_prior: Default False.
    :param independent_filtering: To let DESeq2 perform the independentFiltering or not.
    :param cooks_cutoff: To let DESeq2 decide for the CooksCutoff or specifying a value.
    :param threads: Number of threads to use
    """

    template_file = "rnadiff_light_template.R"
    template_env = Environment(loader=PackageLoader("sequana", "resources/templates"))
    template = template_env.get_template(template_file)

    def __init__(
        self,
        counts_tsv,
        groups_tsv,
        condition,
        comparisons,
        batch,
        fit_type="parametric",
        beta_prior=False,
        independent_filtering=True,
        cooks_cutoff=None,
        threads=4,
    ):
        self.counts_tsv = counts_tsv
        self.groups_tsv = groups_tsv

        self.counts = pd.read_csv(counts_tsv, sep="\t", index_col="Geneid")
        self.groups = pd.read_csv(groups_tsv, sep="\t", index_col="label")

        self.condition = condition
        self.comparisons = comparisons
        self.comparisons_str = (
            f"list({', '.join(['c' + str(x) for x in self.comparisons])})"
        )
        self.batch = batch
        self.design = f"~{batch + '+' + condition if batch else condition}"
        self.fit_type = fit_type
        self.beta_prior = "TRUE" if beta_prior else "FALSE"
        self.independent_filtering = "TRUE" if independent_filtering else "FALSE"
        self.cooks_cutoff = cooks_cutoff if cooks_cutoff else "TRUE"
        self.threads = threads
        self.results = None

    def __repr__(self):
        info = f"RNADiffAnalysis object:\n\
- {self.counts.shape[1]} samples.\n\
- {len(self.comparisons)} comparisons.\n\n\
Counts overview:\n\
{self.counts.head()}\n\n\
Groups overview:\n\
{self.groups.head()}"

        return info

    def run(self):
        """Generate a DESeq2 script from template and execute it."""

        logger.info("Starting differential analysis with DESeq2...")

        with open("rnadiff_light.R", "w") as f:
            f.write(RNADiffAnalysis.template.render(self.__dict__))

        p = subprocess.run(
            ["Rscript", "rnadiff_light.R"], universal_newlines=True, capture_output=True
        )

        # Capture rnadiff output
        with open("rnadiff.err", "w") as f:
            f.write(p.stderr)
        with open("rnadiff.out", "w") as f:
            f.write(p.stdout)

        self.results = [
            pd.read_csv(f"{x}_VS_{y}_degs_DESeq2.csv", index_col=0)
            for x, y in self.comparisons
        ]

        logger.info("DONE")

    def annotate(self, gff):
        """Add annotation to deseq results (from mart or gff)
        **IN DEV**
        """
        gff = GFF3(gff)
        gff_df = gff.get_df()
        gff_df.to_csv("test.csv")


class RNADiffResults:
    """An object representation of results coming from a RNADiff analysis."""

    def __init__(
        self, filename, alpha=0.05, log2_fc=0, pattern="*complete*xls", sep="\t"
    ):
        """.. rubric:: constructor

        :param rnadiff_results: can be a folder for a simple comparison, or an
            output file containing the rseults of a rnadiff analysis
        :param alpha:
        :param out_dir:
        :param log2_fc: the log2 fold change to apply


        Note that if the index or Name/ID columns are made of ensemble ID, they
        may be preceded by the gene: prefix, which is removed
        """
        if isinstance(filename, RNADiffResults):
            self.df = filename.df.copy()
            self.filename = filename.filename
        elif os.path.isfile(filename):
            # must be a 'complete' file from rnadiff analysis
            self.filename = filename
            self.df = pd.read_csv(self.filename, sep, index_col=0)
        elif os.path.isdir(filename):
            filenames = glob.glob(filename + "/" + pattern)
            found_unique = False
            if len(filenames) == 1:
                self.filename = filenames[0]
                found_unique = True
            else:
                filenames = glob.glob(filename + "/tables/" + pattern)
                if len(filenames) == 1:
                    found_unique = True
                    self.filename = filenames[0]
            if found_unique is False:
                raise IOError(
                    "Found several file with the {} pattern. Please be"
                    " more restrictive using the pattern argument"
                )
            self.df = pd.read_csv(self.filename, sep, index_col=0)
        else:
            raise TypeError("{} does not exists".format(filename))

        # Some cleanup
        self.df.index = [x.replace("gene:", "") for x in self.df.index]
        for col in ["ID", "Name"]:
            if col in self.df.columns:
                try:
                    self.df[col] = [x.replace("gene:", "") for x in self.df[col]]
                except:
                    pass

        # Just an alias to a subset of the dataframe
        normed = [x for x in self.df.columns if x.startswith("norm")]
        self.normcounts = self.df[normed]

        # some parameters/attributes
        self._alpha = alpha
        self._log2_fc = log2_fc

        # What are the sample names and condition names
        self.sample_names = [x.replace("norm.", "") for x in normed]
        self.condition_names = set([x[0:-1] for x in self.sample_names])
        self.set_colors()

        self._set_gene_lists_one_condition()

    def set_colors(self, colors=None):
        self.colors = {}
        if colors is None:
            colors = ["orange", "b", "r", "y", "k"]
        for i, name in enumerate(self.condition_names):
            try:
                self.colors[name] = colors[i]
            except:
                self.colors[name] = colors[len(colors) - 1]

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

        gene_sets = {}
        key = "vs".join(sorted(self.condition_names))
        gene_sets = {}

        condition_up = np.logical_and(
            self.df.log2FoldChange > self.log2_fc, self.df.padj < self.alpha
        )
        gene_sets["up"] = set(self.df[condition_up].index)

        condition_down = np.logical_and(
            self.df.log2FoldChange < -self.log2_fc, self.df.padj < self.alpha
        )
        gene_sets["down"] = set(self.df[condition_down].index)

        condition_all = np.logical_or(
            np.logical_and(
                self.df.padj < self.alpha, self.df.log2FoldChange > self.log2_fc
            ),
            np.logical_and(
                self.df.padj < self.alpha, self.df.log2FoldChange < -self.log2_fc
            ),
        )
        gene_sets["all"] = set(self.df[condition_all].index)

        self.gene_lists = gene_sets

    def summary(self):
        """Get a summary DataFrame from a RNADiff analysis."""
        summary = pd.DataFrame(
            {(x, len(self.gene_lists[x])) for x in self.gene_lists.keys()}
        )

        df = summary.set_index(0)
        df.columns = ["_vs_".join(self.condition_names)]
        return df

    def plot_count_per_sample(self, fontsize=12, sample_list=None):
        """ "Number of mapped reads per sample. Each color for each replicate

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

        pylab.bar(
            range(N),
            (dd / 1000000).values,
            color=colors,
            alpha=1,
            zorder=10,
            lw=1,
            ec="k",
            width=0.9,
        )
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

        data = (self.df[self.sample_names] == 0).sum()
        data = data / len(self.df) * 100

        all_null = (self.df[self.sample_names].sum(axis=1) == 0).sum()

        colors = []
        for sample in self.sample_names:
            colors.append(self.colors[self.get_cond_from_sample(sample)])

        pylab.clf()
        pylab.bar(
            range(N), data, color=colors, alpha=1, zorder=10, lw=1, ec="k", width=0.9
        )
        pylab.axhline(
            all_null / len(self.df) * 100, lw=2, ls="--", color="k", zorder=20
        )
        pylab.xticks(range(N), self.sample_names)
        pylab.xlabel("Sample")
        pylab.ylabel("Proportion of null counts (%)")
        pylab.grid(True, zorder=0)

    def plot_volcano(
        self,
        padj=0.05,
        add_broken_axes=False,
        markersize=4,
        limit_broken_line=[20, 40],
        plotly=False,
    ):
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

        if plotly:
            from plotly import express as px

            df = self.df.copy()
            df["log_adj_pvalue"] = -pylab.log10(self.df.padj)
            df["significance"] = ["<0.05" if x else ">=0.05" for x in df.padj < 0.05]

            if "Name" in self.df.columns:
                hover_name = "Name"
            else:
                hover_name = "ID"
            fig = px.scatter(
                df,
                x="log2FoldChange",
                y="log_adj_pvalue",
                hover_name=hover_name,
                log_y=False,
                color="significance",
                height=600,
                labels={"log_adj_pvalue": "log adjusted p-value"},
            )

            return fig

        from brokenaxes import brokenaxes

        M = max(-pylab.log10(self.df.padj.dropna()))

        br1, br2 = limit_broken_line
        if M > br1:
            if add_broken_axes:
                bax = brokenaxes(ylims=((0, br1), (M - 10, M)), xlims=None)
            else:
                bax = pylab
        else:
            bax = pylab

        bax.plot(
            d1.log2FoldChange,
            -np.log10(d1.padj),
            marker="o",
            alpha=0.5,
            color="k",
            lw=0,
            markersize=markersize,
        )
        bax.plot(
            d2.log2FoldChange,
            -np.log10(d2.padj),
            marker="o",
            alpha=0.5,
            color="r",
            lw=0,
            markersize=markersize,
        )

        bax.grid(True)
        try:
            bax.set_xlabel("fold change")
            bax.set_ylabel("log10 adjusted p-value")
        except:
            bax.xlabel("fold change")
            bax.ylabel("log10 adjusted p-value")

        m1 = abs(min(self.df.log2FoldChange))
        m2 = max(self.df.log2FoldChange)
        limit = max(m1, m2)
        try:
            bax.set_xlim([-limit, limit])
        except:
            bax.xlim([-limit, limit])
        try:
            y1, _ = bax.get_ylim()
            ax1 = bax.axs[0].set_ylim([br2, y1[1] * 1.1])
        except:
            y1, y2 = bax.ylim()
            bax.ylim([0, y2])
        bax.axhline(
            -np.log10(0.05), lw=2, ls="--", color="r", label="pvalue threshold (0.05)"
        )
        return bax

    def plot_pca(self, n_components=2, colors=None, plotly=False):
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

        if plotly is True:
            assert n_components == 3
            variance = p.plot(n_components=n_components, colors=colors, show_plot=False)
            from plotly import express as px

            df = pd.DataFrame(p.Xr)
            df.columns = ["PC1", "PC2", "PC3"]
            df["names"] = self.sample_names
            df["colors"] = [colors[x] for x in self.sample_names]
            df["size"] = [10] * len(df)
            df["condition"] = [
                self.get_cond_from_sample(sample) for sample in self.sample_names
            ]
            fig = px.scatter_3d(
                df,
                x="PC1",
                y="PC2",
                z="PC3",
                color="condition",
                labels={
                    "PC1": "PC1 ({}%)".format(round(100 * variance[0], 2)),
                    "PC2": "PC2 ({}%)".format(round(100 * variance[1], 2)),
                    "PC3": "PC3 ({}%)".format(round(100 * variance[2], 2)),
                },
                height=800,
                text="names",
            )
            return fig
        else:
            variance = p.plot(n_components=n_components, colors=colors)

        return variance

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
            import re

            candidates = [
                x
                for x in self.condition_names
                if re.findall("^{}\d+".format(x), sample_name)
            ]
            if len(candidates) == 1:
                return candidates[0]
            else:
                raise ValueError("ambiguous sample name found in several conditions")
        except:
            logger.warning("{} not found or ambiguous name".format(sample_name))
            return None

    def plot_density(self):
        import seaborn

        seaborn.set()
        for x in self.sample_names:
            seaborn.kdeplot(pylab.log10(self.df[x].clip(lower=1)))
            pylab.ylabel("Density")
            pylab.xlabel("Raw counts (log10)")

    def plot_feature_most_present(self):
        """"""
        _description = "test me"
        N = len(self.sample_names)

        data = self.df[self.sample_names].max() / self.df[self.sample_names].sum()
        data = data * 100

        colors = []
        names = []
        for sample in self.sample_names:
            colors.append(self.colors[self.get_cond_from_sample(sample)])

        pylab.clf()
        pylab.barh(
            range(0, len(data)),
            data,
            color=colors,
            alpha=1,
            zorder=10,
            lw=1,
            ec="k",
            height=0.9,
        )
        count = 0
        for x, sample in zip(data, self.sample_names):
            index_name = self.df[sample].argmax()
            try:
                name = self.df.iloc[
                    index_name
                ].Name  # !! this is the index not the column called Name
            except:
                name = self.df.iloc[index_name][
                    "ID"
                ]  # !! this is the index not the column called Name

            pylab.text(
                x, count, name, fontdict={"ha": "center", "va": "center"}, zorder=20
            )
            count += 1
        x0, x1 = pylab.xlim()
        pylab.xlim([0, min(100, x1 * 1.2)])
        pylab.yticks(range(N), self.sample_names)
        pylab.xlabel("Reads captured by most important feature (%s)")
        pylab.ylabel("Sample")

    def plot_dendogram(self, max_features=5000, transform_method="log"):

        assert transform_method in ["log", "anscombe"]
        # first we take the normalised data
        from sequana.viz import clusterisation
        from sequana.viz import dendogram

        cluster = clusterisation.Cluster(self.normcounts)
        # cluster = clusterisation.Cluster(self.df[self.sample_names])
        data = cluster.scale_data(
            transform_method=transform_method, max_features=max_features
        )
        df = pd.DataFrame(data[0])
        df.index = data[1]
        df.columns = self.normcounts.columns

        d = dendogram.Dendogram(df.T)
        d.category = {}
        conditions = list(self.condition_names)
        for sample_name in df.columns:
            cond = self.get_cond_from_sample(sample_name.replace("norm.", ""))
            d.category[sample_name] = conditions.index(cond)
        d.plot()

    def boxplot_rawdata(self, fliersize=2, linewidth=2, **kwargs):
        import seaborn as sbn

        sbn.set()

        colors = []
        for sample in self.sample_names:
            col = self.colors[self.get_cond_from_sample(sample)]
            if col == "orange":
                colors.append("#FFA500")
            else:
                colors.append("#3c5688")

        ax = sbn.boxplot(
            data=self.df[self.sample_names].clip(1),
            linewidth=linewidth,
            fliersize=fliersize,
            palette=colors,
            **kwargs,
        )
        ax.set(yscale="log")

    def boxplot_normeddata(self, fliersize=2, linewidth=2, **kwargs):
        import seaborn as sbn

        sbn.set()
        colors = []
        for sample in self.sample_names:
            col = self.colors[self.get_cond_from_sample(sample)]
            if col == "orange":
                colors.append("#FFA500")
            else:
                colors.append("#3c5688")
        ax = sbn.boxplot(
            data=self.normcounts.clip(1),
            linewidth=linewidth,
            fliersize=fliersize,
            palette=colors,
            **kwargs,
        )
        ax.set(yscale="log")

    def plot_pvalue_hist(self, bins=60, fontsize=16, rotation=0):
        pylab.hist(self.df.pvalue.dropna(), bins=bins, ec="k")
        pylab.grid(True)
        pylab.xlabel("raw p-value", fontsize=fontsize)
        pylab.ylabel("Occurences", fontsize=fontsize)
        try:
            pylab.tight_layout()
        except:
            pass

    def plot_padj_hist(self, bins=60, fontsize=16):
        pylab.hist(self.df.padj.dropna(), bins=bins, ec="k")
        pylab.grid(True)
        pylab.xlabel("Adjusted p-value", fontsize=fontsize)
        pylab.ylabel("Occurences", fontsize=fontsize)
        try:
            pylab.tight_layout()
        except:
            pass

    def plot_dispersion(self):

        pylab.plot(
            self.normcounts.mean(axis=1),
            self.df.dispGeneEst,
            "ok",
            label="Estimate",
            ms=1,
        )
        pylab.plot(
            self.normcounts.mean(axis=1), self.df.dispersion, "ob", label="final", ms=1
        )
        pylab.plot(
            self.normcounts.mean(axis=1), self.df.dispFit, "or", label="Fit", ms=1
        )
        pylab.legend()
        ax = pylab.gca()
        ax.set(yscale="log")
        ax.set(xscale="log")
        pylab.ylabel("Dispersion", fontsize=fontsize)
        pylab.xlabel("Mean of normalized counts", fontsize=fontsize)
