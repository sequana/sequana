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
from jinja2 import Environment, PackageLoader
import subprocess
import seaborn as sns
import matplotlib.pyplot as plt
from itertools import combinations
import os

from sequana.lazy import pandas as pd
from sequana.lazy import pylab
from sequana.lazy import numpy as np

from sequana.gff3 import GFF3
from sequana.viz import Volcano
from sequana.enrichment import PantherEnrichment, KeggPathwayEnrichment

import colorlog
logger = colorlog.getLogger(__name__)


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
    :param gff: Path to the corresponding gff3 to add annotations.
    :param fc_attribute: GFF attribute used in FeatureCounts.
    :param fc_feature: GFF feaure used in FeatureCounts.
    :param annot_cols: GFF attributes to use for results annotations
    :param threads: Number of threads to use
    :param outdir: Path to output directory.
    :param sep: The separator to use in dataframe exports
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
        gff=None,
        fc_attribute=None,
        fc_feature=None,
        annot_cols=["ID", "Name", "gene_biotype"],
        threads=4,
        outdir="rnadiff",
        sep_counts="\t",
        sep_groups="\t",
    ):

        self.counts_tsv = counts_tsv
        self.groups_tsv = groups_tsv

        self.counts = pd.read_csv(counts_tsv, sep=sep_counts, index_col="Geneid",
            comment="#")
        self.groups = pd.read_csv(groups_tsv, sep=sep_groups, index_col="label",
            comment="#")


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
        self.gff = gff
        self.fc_feature = fc_feature
        self.fc_attribute = fc_attribute
        self.threads = threads

        self.outdir = Path(outdir)
        self.results = self._run()

    def __repr__(self):
        info = f"RNADiffAnalysis object:\n\
- {self.counts.shape[1]} samples.\n\
- {len(self.comparisons)} comparisons.\n\n\
Counts overview:\n\
{self.counts.head()}\n\n\
Groups overview:\n\
{self.groups.head()}"

        return info

    def _run(self):
        """Create outdir and a DESeq2 script from template for analysis. Then execute
        this script.
        """

        try:
            self.outdir.mkdir()
        except:
            pass

        with open("rnadiff_light.R", "w") as f:
            f.write(RNADiffAnalysis.template.render(self.__dict__))

        logger.info("Starting differential analysis with DESeq2...")

        p = subprocess.run(
            ["Rscript", "rnadiff_light.R"], universal_newlines=True, capture_output=True
        )

        # Capture rnadiff output
        with open("rnadiff.err", "w") as f:
            f.write(p.stderr)
        with open("rnadiff.out", "w") as f:
            f.write(p.stdout)

        return RNADiffResults(
            self.outdir,
            self.groups_tsv,
            group=self.condition,
            gff=self.gff,
            fc_feature=self.fc_feature,
            fc_attribute=self.fc_attribute,
        )


class RNADiffTable:
    def __init__(self, path, gff=None, alpha=0.05, log2_fc=0, sep="\t"):
        """ A representation of the results of a single rnadiff comparison """
        self.path = Path(path)
        self.name = self.path.stem.replace("_degs_DESeq2", "").replace("-", "_")

        self._alpha = alpha
        self._log2_fc = log2_fc

        self.gff = gff

        self.df = pd.read_csv(self.path, index_col=0, sep=sep)

        self.filt_df = self.filter()
        self.set_gene_lists()

    @property
    def alpha(self):
        return self._alpha

    @alpha.setter
    def alpha(self, value):
        self._alpha = value
        self.filt_df = self.filter()
        self.set_gene_lists()

    @property
    def log2_fc(self):
        return self._log2_fc

    @log2_fc.setter
    def log2_fc(self, value):
        self._log2_fc = value
        self.filt_df = self.filter()
        self.set_gene_lists()

    def filter(self):
        """filter a DESeq2 result with FDR and logFC thresholds"""

        fc_filt = self.df["log2FoldChange"].abs() >= self._log2_fc
        fdr_filt = self.df["padj"] <= self._alpha

        return self.df[fc_filt.values & fdr_filt.values]

    def set_gene_lists(self):

        self.gene_lists = {
            "up": list(self.filt_df.query("log2FoldChange > 0").index),
            "down": list(self.filt_df.query("log2FoldChange < 0").index),
            "all": list(self.filt_df.index),
        }

    def summary(self):

        return pd.DataFrame(
            {
                "log2_fc": self._log2_fc,
                "alpha": self._alpha,
                "up": len(self.gene_lists["up"]),
                "down": len(self.gene_lists["down"]),
                "all": len(self.gene_lists["all"]),
            },
            index=[self.name],
        )

    def plot_volcano(self):

        fc = self.df.loc[:, "log2FoldChange"]
        pvalues = self.df.loc[:, "padj"]
        Volcano(fc, pvalues).plot()
        plt.title(self.name)

    def plot_pvalue_hist(self, bins=60, fontsize=16, rotation=0):

        plt.hist(self.df.pvalue.dropna(), bins=bins, ec="k")
        plt.xlabel("raw p-value")
        plt.ylabel("Occurences")


class RNADiffResults:
    def __init__(
        self,
        rnadiff_folder,
        meta,
        gff=None,
        fc_attribute=None,
        fc_feature=None,
        pattern="*vs*_degs_DESeq2.tsv",
        alpha=0.05,
        log2_fc=0,
        sep="\t",
        palette=sns.color_palette(desat=0.6),
        group="condition",
        annot_cols=["ID", "Name", "gene_biotype"],
    ):
        """
        :rnadiff_folder:
        """
        self.path = Path(rnadiff_folder)
        self.files = [x for x in self.path.glob(pattern)]

        self.meta = self._get_meta(meta, sep=sep, group=group, palette=palette)

        self.counts_raw = pd.read_csv(
            self.path / "counts_raw.tsv", index_col=0, sep=sep
        )
        self.counts_norm = pd.read_csv(
            self.path / "counts_normed.tsv", index_col=0, sep=sep
        )
        self.counts_vst = pd.read_csv(
            self.path / "counts_vst_norm.tsv", index_col=0, sep=sep
        )
        self.dds_stats = pd.read_csv(
            self.path / "overall_dds.tsv", index_col=0, sep=sep
        )

        self.gff = gff
        self.fc_attribute = fc_attribute
        self.fc_feature = fc_feature
        self.annot_cols = annot_cols

        self._alpha = alpha
        self._log2_fc = log2_fc

        if self.gff:
            self.annot_df = self._get_annot()

        self.comparisons = self.import_tables()

        self.df = self._get_total_df()
        self.filt_df = self._get_total_df(filtered=True)

    @property
    def alpha(self):
        return self._alpha

    @alpha.setter
    def alpha(self, value):
        self._alpha = value
        self.comparisons = self.import_tables()
        self.filt_df = self._get_total_df(filtered=True)

    @property
    def log2_fc(self):
        return self._log2_fc

    @log2_fc.setter
    def log2_fc(self, value):
        self._log2_fc = value
        self.comparisons = self.import_tables()
        self.filt_df = self._get_total_df(filtered=True)

    def import_tables(self):

        return {
            compa.stem.replace("_degs_DESeq2", "").replace("-", "_"): RNADiffTable(
                compa, alpha=self._alpha, log2_fc=self._log2_fc
            )
            for compa in self.files
        }

    def _get_annot(
        self,
    ):
        """Get a properly formatted dataframe from the gff."""

        df = GFF3(self.gff).get_df()
        df = df.query("type == @self.fc_feature").loc[:, self.annot_cols]
        df.drop_duplicates(inplace=True)
        df.set_index(self.fc_attribute, inplace=True)
        df.columns = pd.MultiIndex.from_product([["annotation"], df.columns])

        return df

    def _get_total_df(self, filtered=False):
        """Concatenate all rnadiff results in a single dataframe.

        FIXME: Columns relative to significative comparisons are not using
        self.log2_fc and self.alpha
        """

        dfs = []

        for compa, res in self.comparisons.items():
            df = res.filt_df if filtered else res.df
            df = df.transpose().reset_index()
            df["file"] = res.name
            df = df.set_index(["file", "index"])
            dfs.append(df)

        df = pd.concat(dfs, sort=True).transpose()

        # Add number of comparisons which are significative for a given gene
        num_sign_compa = (df.loc[:, (slice(None), "padj")] < 0.05).sum(axis=1)
        df.loc[:, ("statistics", "num_of_significative_comparisons")] = num_sign_compa

        # Add list of comparisons which are significative for a given gene
        df_sign_padj = df.loc[:, (slice(None), "padj")] < 0.05
        sign_compa = df_sign_padj.loc[:, (slice(None), "padj")].apply(
            # Extract column names (comparison names) for significative comparisons
            lambda row: {col_name[0] for sign, col_name in zip(row, row.index) if sign},
            axis=1,
        )
        df.loc[:, ("statistics", "significative_comparisons")] = sign_compa

        if self.gff and self.fc_attribute and self.fc_feature:
            df = pd.concat([self.annot_df, df], axis=1)
        else:
            logger.warning(
                "Missing any of gff, fc_attribute or fc_feature. No annotation will be added."
            )

        return df

    def summary(self):
        return pd.concat(res.summary() for compa, res in self.comparisons.items())

    def report(self):

        template_file = "rnadiff_report.html"
        template_env = Environment(
            loader=PackageLoader("sequana", "resources/templates")
        )
        template = template_env.get_template(template_file)

        with open("rnadiff_report.html", "w") as f:
            f.write(
                template.render(
                    {"table": self.summary().to_html(classes="table table-striped")}
                )
            )

    def run_enrichment_go(self, taxon, annot_col="Name", out_dir="enrichment"):

        out_dir = Path(out_dir) / "figures"
        out_dir.mkdir(exist_ok=True, parents=True)

        gene_lists_dict = self.get_gene_lists(annot_col=annot_col, Nmax=2000)
        enrichment = {}
        ontologies = {"GO:0003674": "BP", "GO:0008150": "MF", "GO:0005575": "CC"}

        for compa in self.comparisons:
            gene_lists = gene_lists_dict[compa]
            pe = PantherEnrichment(gene_lists, taxon)
            pe.compute_enrichment(ontologies=ontologies.keys(), progress=False)

            for direction in ["up", "down", "all"]:
                for ontology in ontologies.keys():
                    enrichment[(compa, direction, ontology)] = pe.get_data(
                        direction, ontology, include_negative_enrichment=False
                    )
                    plt.figure()
                    try:
                        pe.plot_go_terms(direction, ontology, compute_levels=False)
                    except:
                        return pe
                    plt.tight_layout()
                    plt.savefig(
                        out_dir / f"go_{compa}_{direction}_{ontologies[ontology]}.pdf"
                    )

            logger.info(f"Panther enrichment for {compa} DONE.")

        df = pd.concat(enrichment).sort_index()
        df.index.rename(
            ["comparison", "direction", "GO_category", "index"], inplace=True
        )

        self.enrichment_go = df

        # Export results (should be moved to enrichment.py at some point I think)
        with pd.ExcelWriter(out_dir.parent / "enrichment_go.xlsx") as writer:
            df = self.enrichment_go.copy()
            df.reset_index(inplace=True)
            df.to_excel(writer, "go", index=False)
            ws = writer.sheets["go"]
            ws.autofilter(0, 0, df.shape[0], df.shape[1] - 1)

    def run_enrichment_kegg(self, organism, annot_col="Name", out_dir="enrichment"):

        out_dir = Path(out_dir) / "figures"
        out_dir.mkdir(exist_ok=True, parents=True)

        gene_lists_dict = self.get_gene_lists(annot_col=annot_col)
        enrichment = {}

        for compa in self.comparisons:
            gene_lists = gene_lists_dict[compa]
            ke = KeggPathwayEnrichment(gene_lists, organism, progress=False)
            ke.compute_enrichment()

            for direction in ["up", "down", "all"]:
                enrichment[(compa, direction)] = ke._get_final_df(
                    ke.enrichment[direction].results, nmax=10000
                )
                plt.figure()
                ke.scatterplot(direction)
                plt.tight_layout()
                plt.savefig(out_dir / f"kegg_{compa}_{direction}.pdf")

            logger.info(f"KEGG enrichment for {compa} DONE.")

        df = pd.concat(enrichment).sort_index()
        df.index.rename(["comparison", "direction", "index"], inplace=True)

        self.enrichment_kegg = df

        # Export results (should be moved to enrichment.py at some point I think)
        with pd.ExcelWriter(out_dir.parent / "enrichment_kegg.xlsx") as writer:
            df = self.enrichment_kegg.copy()
            df.reset_index(inplace=True)
            df.to_excel(writer, "kegg", index=False)
            ws = writer.sheets["kegg"]
            ws.autofilter(0, 0, df.shape[0], df.shape[1] - 1)

    def get_gene_lists(self, annot_col="index", Nmax=None):

        gene_lists_dict = {}

        for compa in self.comparisons.keys():
            df = self.df.loc[:, ["annotation", compa]].copy()
            df = df.droplevel(0, axis=1)

            fc_filt = df["log2FoldChange"].abs() >= self._log2_fc
            fdr_filt = df["padj"] <= self._alpha

            df = df[fc_filt.values & fdr_filt.values]
            df.reset_index(inplace=True)

            if Nmax:
                df.sort_values("log2FoldChange", ascending=False, inplace=True)
                up_genes = list(df.query("log2FoldChange > 0")[annot_col])[:Nmax]

                df.sort_values("log2FoldChange", ascending=True, inplace=True)
                down_genes = list(df.query("log2FoldChange < 0")[annot_col])[:Nmax]

                all_genes = list(
                    list(
                        df.sort_values("log2FoldChange", key=abs, ascending=False)[
                            annot_col
                        ]
                    )[:Nmax]
                )

            else:
                up_genes = list(df.query("log2FoldChange > 0")[annot_col])
                down_genes = list(df.query("log2FoldChange < 0")[annot_col])
                all_genes = list(df.loc[:, annot_col])

            gene_lists_dict[compa] = {
                "up": up_genes,
                "down": down_genes,
                "all": all_genes,
            }

        return gene_lists_dict

    def _get_meta(self, meta_file, sep, group, palette):
        """Import metadata from a table file and add color groups following the
        groups defined in the column 'group' of the table file.
        """
        meta_df = pd.read_csv(meta_file, sep=sep, index_col=0)
        col_map = dict(zip(meta_df.loc[:, group].unique(), palette))
        meta_df["group_color"] = meta_df.loc[:, group].map(col_map)

        return meta_df

    def _format_plot(self, title, xlabel, ylabel):
        plt.title(title)
        plt.xticks(rotation=45, ha="right")
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)

    def get_specific_commons(self, direction, compas=None, annot_col="index"):
        """From all the comparisons contained by the object, extract gene lists which
        are common (but specific, ie a gene appear only appears in the
        combination considered) comparing all combinations of comparisons.

        :param direction: The regulation direction (up, down or all) of the gene
        lists to consider

        :param compas: Specify a list of comparisons to consider (Comparisons
        names can be found with self.comparisons.keys()).

        """

        common_specific_dict = {}

        total_gene_lists = self.get_gene_lists(annot_col=annot_col)

        if not compas:
            compas = self.comparisons.keys()

        for size in range(1, len(compas)):
            for compa_group in combinations(compas, size):
                gene_lists = [
                    total_gene_lists[compa][direction] for compa in compa_group
                ]
                commons = set.intersection(
                    *[set(gene_list) for gene_list in gene_lists]
                )
                other_compas = [compa for compa in compas if compa not in compa_group]
                genes_in_other_compas = {
                    x
                    for other_compa in other_compas
                    for x in total_gene_lists[other_compa][direction]
                }

                commons = commons - genes_in_other_compas
                common_specific_dict[compa_group] = commons

        return common_specific_dict

    def plot_count_per_sample(self):
        """Number of mapped and annotated reads (i.e. counts) per sample. Each color
        for each replicate
        """

        df = self.counts_raw.sum().rename("total_counts")
        df = pd.concat([self.meta, df], axis=1)

        plt.bar(df.index, df.total_counts, color=df.group_color)

        self._format_plot(
            title="Total counts", xlabel="Sample", ylabel="Total number of counts"
        )

    def plot_percentage_null_read_counts(self):
        """Bars represent the percentage of null counts in each samples.  The dashed
        horizontal line represents the percentage of feature counts being equal
        to zero across all samples"""

        df = (self.counts_raw == 0).sum() / self.counts_raw.shape[0] * 100
        df = df.rename("percent_null")
        df = pd.concat([self.meta, df], axis=1)

        plt.bar(df.index, df.percent_null, color=df.group_color)

        all_null = (self.counts_raw == 0).all(axis=1).sum() / self.counts_raw.shape[0]

        plt.axhline(all_null, ls="--", color="black", alpha=0.5)

        plt.xticks(rotation=45, ha="right")
        plt.xlabel("Sample")
        plt.ylabel("Proportion of null counts (%)")

        self._format_plot(
            title="Proportion of null counts",
            xlabel="Sample",
            ylabel="% of null counts",
        )

    def plot_pca(
        self,
        n_components=2,
        colors=None,
        plotly=False,
        max_features=500,
        genes_to_remove=[],
    ):

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

        # Get most variable genes (n=max_features)
        top_features = (
            self.counts_vst.var(axis=1)
            .sort_values(ascending=False)
            .index[:max_features]
        )

        if genes_to_remove:
            top_features = [x for x in top_features if x not in genes_to_remove]

        counts_top_features = self.counts_vst.loc[top_features, :]

        p = PCA(counts_top_features)

        if plotly is True:
            assert n_components == 3
            variance = p.plot(
                n_components=n_components,
                colors=colors,
                show_plot=False,
                max_features=max_features,
            )
            from plotly import express as px

            df = pd.DataFrame(p.Xr)
            df.columns = ["PC1", "PC2", "PC3"]
            df["size"] = [10] * len(df)
            df = pd.concat([df, self.meta], axis=1, ignore_index=True)
            return df

            fig = px.scatter_3d(
                df,
                x="PC1",
                y="PC2",
                z="PC3",
                color="group_color",
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
            variance = p.plot(
                n_components=n_components,
                colors=self.meta.group_color,
                max_features=max_features,
            )

        return variance

    def plot_mds(self, n_components=2, colors=None, clf=True):
        """IN DEV, not functional"""

        from sequana.viz.mds import MDS

        p = MDS(self.counts[self.sample_names])
        if colors is None:
            colors = {}
            for sample in self.sample_names:
                colors[sample] = self.colors[self.get_cond_from_sample(sample)]
        p.plot(n_components=n_components, colors=colors, clf=clf)

    def plot_isomap(self, n_components=2, colors=None):
        """IN DEV, not functional"""

        from sequana.viz.isomap import Isomap

        p = Isomap(self.df[self.sample_names])
        if colors is None:
            colors = {}
            for sample in self.sample_names:
                colors[sample] = self.colors[self.get_cond_from_sample(sample)]
        p.plot(n_components=n_components, colors=colors)

    def plot_density(self):
        import seaborn

        seaborn.set()
        for sample in self.counts_raw.columns:
            seaborn.kdeplot(pylab.log10(self.counts_raw[sample].clip(lower=1)))

        self._format_plot(
            title="Count density distribution",
            xlabel="Raw counts (log10)",
            ylabel="Density",
        )

    def plot_feature_most_present(self):
        """"""
        description = "test me"

        df = []

        for x, y in self.counts_raw.idxmax().iteritems():

            most_exp_gene_count = self.counts_raw.stack().loc[y, x]
            total_sample_count = self.counts_raw.sum().loc[x]

            df.append(
                {
                    "label": x,
                    "gene_id": y,
                    "count": most_exp_gene_count,
                    "total_sample_count": total_sample_count,
                    "most_exp_percent": most_exp_gene_count / total_sample_count * 100,
                }
            )

        df = pd.DataFrame(df).set_index("label")
        df = pd.concat([self.meta, df], axis=1)

        p = plt.bar(df.index, df.most_exp_percent, color=df.group_color)

        for idx, rect in enumerate(p):
            plt.text(
                rect.get_x() + rect.get_width() / 2.0,
                0.95 * rect.get_height(),
                df.gene_id.iloc[idx],
                ha="center",
                va="top",
                rotation=90,
            )

        self._format_plot(
            title="Counts monopolized by the most expressed gene",
            xlabel="Sample",
            ylabel="Percent of total reads",
        )

    def plot_dendogram(
        self,
        max_features=5000,
        transform_method="log",
        method="ward",
        metric="euclidean",
    ):
        # for info about metric and methods: https://tinyurl.com/yyhk9cl8

        assert transform_method in ["log", "anscombe", None]
        # first we take the normalised data
        from sequana.viz import clusterisation
        from sequana.viz import dendogram

        cluster = clusterisation.Cluster(self.counts_norm)
        # cluster = clusterisation.Cluster(self.df[self.sample_names])
        if transform_method is not None:
            data = cluster.scale_data(
                transform_method=transform_method, max_features=max_features
            )
            df = pd.DataFrame(data[0])
            df.index = data[1]
            df.columns = self.counts_norm.columns
        else:
            df = pd.DataFrame(self.counts_norm)
            # df.index = data[1]
            df.columns = self.counts_norm.columns

        d = dendogram.Dendogram(
            df.T,
            metric=metric,
            method=method,
            side_colors=list(self.meta.group_color.unique()),
        )

        # Convert groups into numbers for Dendrogram category
        group_conv = {group: i for i, group in enumerate(self.meta.condition.unique())}
        d.category = self.meta.condition.map(group_conv).to_dict()
        d.plot()

    def plot_boxplot_rawdata(self, fliersize=2, linewidth=2, **kwargs):
        import seaborn as sbn

        ax = sbn.boxplot(
            data=self.counts_raw.clip(1),
            linewidth=linewidth,
            fliersize=fliersize,
            palette=self.meta.group_color,
            **kwargs,
        )
        ax.set_yscale("log")
        self._format_plot(
            title="Raw count distribution", xlabel="Samples", ylabel="Raw counts"
        )

    def plot_boxplot_normeddata(self, fliersize=2, linewidth=2, **kwargs):
        import seaborn as sbn

        ax = sbn.boxplot(
            data=self.counts_norm.clip(1),
            linewidth=linewidth,
            fliersize=fliersize,
            palette=self.meta.group_color,
            **kwargs,
        )
        ax.set(yscale="log")

    def plot_dispersion(self):

        pylab.plot(
            self.dds_stats.baseMean,
            self.dds_stats.dispGeneEst,
            "ok",
            label="Estimate",
            ms=1,
        )
        pylab.plot(
            self.dds_stats.baseMean,
            self.dds_stats.dispersion,
            "ob",
            label="final",
            ms=1,
        )
        pylab.plot(
            self.dds_stats.baseMean, self.dds_stats.dispFit, "or", label="Fit", ms=1
        )
        pylab.legend()
        ax = pylab.gca()
        ax.set(yscale="log")
        ax.set(xscale="log")

        self._format_plot(
            title="Dispersion estimation",
            xlabel="Mean of normalized counts",
            ylabel="Dispersion",
        )
