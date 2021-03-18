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
import sys
import os

from pathlib import Path
from jinja2 import Environment, PackageLoader
import subprocess
import seaborn as sns
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

__all__ = ["RNADiffAnalysis", "RNADiffResults", "RNADiffTable", "RNADesign"]


class RNADesign:
    """Simple RNA design handler"""

    def __init__(self, filename, reference=None):
        self.filename = filename
        self.df = pd.read_csv(filename, sep=",")
        self.reference = reference

    def _get_conditions(self):
        return sorted(self.df.condition.unique())

    conditions = property(_get_conditions)

    def _get_comparisons(self):
        conditions = self.conditions
        if self.reference is None:
            import itertools

            comps = list(itertools.combinations(conditions, 2))
        else:
            # only those versus reference
            comps = [(x, self.reference) for x in conditions if x != self.reference]
        return sorted(comps)

    comparisons = property(_get_comparisons)


class RNADiffAnalysis:
    """A tool to prepare and run a RNA-seq differential analysis with DESeq2

    :param counts_file: Path to tsv file out of FeatureCount with all samples together.
    :param design_file: Path to tsv file with the definition of the groups for each sample.
    :param condition: The name of the column from groups_tsv to use as condition. For more
        advanced design, a R function of the type 'condition*inter' (without the '~') could
        be specified (not tested yet). Each name in this function should refer to column
        names in groups_tsv.
    :param comparisons: A list of tuples indicating comparisons to be made e.g A vs B would be [("A", "B")]
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

    This class reads a :class:`sequana.featurecounts.`

    r = rnadiff.RNADiffAnalysis("counts.csv", "design.csv",
            condition="condition", comparisons=[(("A", "B"), ('A', "C")],
            fc_feature="gene",
            fc_attribute="ID", gff="mygff.gff")


    """

    _template_file = "rnadiff_light_template.R"
    _template_env = Environment(loader=PackageLoader("sequana", "resources/templates"))
    template = _template_env.get_template(_template_file)

    def __init__(
        self,
        counts_file,
        design_file,
        condition,
        comparisons,
        batch=None,
        fit_type="parametric",
        beta_prior=False,
        independent_filtering=True,
        cooks_cutoff=None,
        gff=None,
        fc_attribute=None,
        fc_feature=None,
        annot_cols=None,
        # annot_cols=["ID", "Name", "gene_biotype"],
        threads=4,
        outdir="rnadiff",
        sep_counts=",",
        sep_design=",",
    ):

        self.counts_filename = counts_file
        self.design_filename = design_file

        self.counts = pd.read_csv(
            counts_file, sep=sep_counts, index_col="Geneid", comment="#"
        )
        self.design = pd.read_csv(
            design_file,
            sep=sep_design,
            comment="#",
            dtype={"label": str},
        ).set_index("label")

        # TODO: Check and resorting if necessary
        if sorted(list(self.counts.columns)) != sorted(list(self.design.index)):
            logger.error(f"Counts columns and design rows does not match.")
            logger.error(self.counts.columns)
            logger.error(self.design.index)
            sys.exit(1)

        # set condition of the statistical model
        columns = ",".join(self.design.columns)
        if condition not in columns:
            logger.error(
                f"""Your condition named '{condition}' is expected to
be in the header of your design file but was not found. Candidates are:
    {columns}"""
            )
            sys.exit(1)
        self.condition = condition
        self.comparisons = comparisons

        # let us check the consistenct of the design and comparisons
        valid_conditions = ",".join(set(self.design[condition].values))
        for item in [x for y in self.comparisons for x in y]:
            if item not in self.design[condition].values:
                logger.error(
                    f"""{item} not found in the design. Fix the design
or comparisons. possible values are {valid_conditions}"""
                )
                sys.exit(1)

        self.comparisons_str = (
            f"list({', '.join(['c' + str(x) for x in self.comparisons])})"
        )
        self.batch = batch
        self.model = f"~{batch + '+' + condition if batch else condition}"
        self.fit_type = fit_type
        self.beta_prior = "TRUE" if beta_prior else "FALSE"
        self.independent_filtering = "TRUE" if independent_filtering else "FALSE"
        self.cooks_cutoff = cooks_cutoff if cooks_cutoff else "TRUE"
        self.gff = gff
        self.fc_feature = fc_feature
        self.fc_attribute = fc_attribute
        self.threads = threads

        self.outdir = Path(outdir)

    def __repr__(self):
        info = f"RNADiffAnalysis object:\n\
- {self.counts.shape[1]} samples.\n\
- {len(self.comparisons)} comparisons.\n\n\
Counts overview:\n\
{self.counts.head()}\n\n\
Design overview:\n\
{self.design.head()}"

        return info

    def run(self):
        """Create outdir and a DESeq2 script from template for analysis. Then execute
        this script.

        :return: a :class:`RNADiffResults` instance
        """
        logger.info("Running DESeq2 analysis. Please wait")
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

        results = RNADiffResults(
            self.outdir,
            self.design_filename,
            group=self.condition,
            gff=self.gff,
            fc_feature=self.fc_feature,
            fc_attribute=self.fc_attribute,
        )
        return results


class RNADiffTable:
    def __init__(self, path, alpha=0.05, log2_fc=0, sep=","):
        """ A representation of the results of a single rnadiff comparison """
        self.path = Path(path)
        self.name = self.path.stem.replace("_degs_DESeq2", "").replace("-", "_")

        self._alpha = alpha
        self._log2_fc = log2_fc

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
            df["significance"] = [
                "<{}".format(padj) if x else ">={}".format(padj) for x in df.padj < padj
            ]

            if "Name" in self.df.columns:
                hover_name = "Name"
            elif "gene_id" in self.df.columns:
                hover_name = "gene_id"
            elif "locus_tag" in self.df.columns:
                hover_name = "locus_tag"
            elif "ID" in self.df.columns:
                hover_name = "ID"
            else:
                hover_name = None
            fig = px.scatter(
                df,
                x="log2FoldChange",
                y="log_adj_pvalue",
                hover_name=hover_name,
                hover_data=["baseMean"],
                log_y=False,
                opacity=0.5,
                color="significance",
                height=600,
                labels={"log_adj_pvalue": "log adjusted p-value"},
            )
            # axes[0].axhline(
            # -np.log10(0.05), lw=2, ls="--", color="r", label="pvalue threshold (0.05)"
            # i)
            # in future version of plotly, a add_hlines will be available. For
            # now, this is the only way to add axhline
            fig.update_layout(
                shapes=[
                    dict(
                        type="line",
                        xref="x",
                        x0=df.log2FoldChange.min(),
                        x1=df.log2FoldChange.max(),
                        yref="y",
                        y0=-pylab.log10(padj),
                        y1=-pylab.log10(padj),
                        line=dict(color="black", width=1, dash="dash"),
                    )
                ]
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

        if colors is None:
            colors = {}
            for sample in self.sample_names:
                colors[sample] = self.colors[self.get_cond_from_sample(sample)]

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
            variance = p.plot(
                n_components=n_components, colors=colors, max_features=max_features
            )

        return variance

    def plot_pvalue_hist(self, bins=60, fontsize=16, rotation=0):

        pylab.hist(self.df.pvalue.dropna(), bins=bins, ec="k")
        pylab.xlabel("raw p-value")
        pylab.ylabel("Occurences")

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


class RNADiffResults:
    """The output of a RNADiff analysis"""

    def __init__(
        self,
        rnadiff_folder,
        design_file=None,
        gff=None,
        fc_attribute=None,
        fc_feature=None,
        pattern="*vs*_degs_DESeq2.csv",
        alpha=0.05,
        log2_fc=0,
        palette=sns.color_palette(desat=0.6),
        group="condition",
        annot_cols=None,
        # annot_cols=["ID", "Name", "gene_biotype"],
    ):
        """

        :rnadiff_folder:

        """
        self.path = Path(rnadiff_folder)
        self.files = [x for x in self.path.glob(pattern)]

        self.counts_raw = pd.read_csv(
            self.path / "counts_raw.csv", index_col=0, sep=","
        )
        self.counts_raw.sort_index(axis=1, inplace=True)

        self.counts_norm = pd.read_csv(
            self.path / "counts_normed.csv", index_col=0, sep=","
        )
        self.counts_norm.sort_index(axis=1, inplace=True)

        self.counts_vst = pd.read_csv(
            self.path / "counts_vst_norm.csv", index_col=0, sep=","
        )
        self.counts_vst.sort_index(axis=1, inplace=True)

        self.dds_stats = pd.read_csv(
            self.path / "overall_dds.csv", index_col=0, sep=","
        )

        # read different results and sort by sample name all inputs
        # TODO make this a function to be reused in RNADiffAnalysis for example
        if design_file == None:
            conditions = []
            labels = self.counts_raw.columns
            for label in labels:
                condition = input(
                    f"Please give use a condition name for the {label} label: "
                )
                conditions.append(condition)
            df = pd.DataFrame({"label": labels, "condition": conditions})
            df.set_index("label", inplace=True)
            df.sort_index(inplace=True)
            col_map = dict(zip(df.loc[:, group].unique(), palette))
            df["group_color"] = df.loc[:, group].map(col_map)
            self.design_df = df
        else:
            self.design_df = self._get_design(design_file, group=group, palette=palette)
            self.design_df.sort_index(inplace=True)
            self.design = RNADesign(design_file)

        # optional annotation
        self.fc_attribute = fc_attribute
        self.fc_feature = fc_feature
        self.annot_cols = annot_cols
        if gff:
            if fc_feature is None or fc_attribute is None:
                logger.error(
                    "Since you provided a GFF filem you must provide the feature and attribute to be used."
                )
            self.annotation = self.read_annot(gff)
        else:
            self.annotation = None

        # some filtering attributes
        self._alpha = alpha
        self._log2_fc = log2_fc

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

    def to_csv(self, filename):
        self.df.to_csv(filename)

    def read_csv(self, filename):
        self.df = pd.read_csv(filename, index_col=0, header=[0, 1])

    def import_tables(self):
        from easydev import AttrDict

        data = {
            compa.stem.replace("_degs_DESeq2", "").replace("-", "_"): RNADiffTable(
                compa,
                alpha=self._alpha,
                log2_fc=self._log2_fc,
                # gff=self.annotation.annotation,
            )
            for compa in self.files
        }

        return AttrDict(**data)

    def read_annot(self, gff_filename):
        """Get a properly formatted dataframe from the gff."""

        gff = GFF3(gff_filename)
        df = gff.get_df()

        if self.annot_cols is None:
            lol = [
                list(x.keys())
                for x in df.query("type==@self.fc_feature")["attributes"].values
            ]
            annot_cols = list(set([x for item in lol for x in item]))
        else:
            annot_cols = self.annot_cols

        df = df.query("type == @self.fc_feature").loc[:, annot_cols]
        df.drop_duplicates(inplace=True)

        df.set_index(self.fc_attribute, inplace=True)

        # It may happen that a GFF has duplicated IDs ! For instance ecoli
        # has 20 duplicated ID that are part 1 and 2 of the same gene
        df = df[~df.index.duplicated(keep="last")]

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

        if self.annotation is not None and self.fc_attribute and self.fc_feature:
            df = pd.concat([self.annotation, df], axis=1)
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

    def run_enrichment_go(
        self, taxon, annot_col="Name", out_dir="enrichment"
    ):  # pragma: no cover

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
                    pylab.figure()
                    try:
                        pe.plot_go_terms(direction, ontology, compute_levels=False)
                    except:
                        return pe
                    pylab.tight_layout()
                    pylab.savefig(
                        out_dir / f"go_{compa}_{direction}_{ontologies[ontology]}.pdf"
                    )
                    pylab.savefig(
                        out_dir / f"go_{compa}_{direction}_{ontologies[ontology]}.png"
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
            try:
                ws.autofilter(0, 0, df.shape[0], df.shape[1] - 1)
            except:
                logger.warning("XLS formatting issue.")

    def run_enrichment_kegg(
        self, organism, annot_col="Name", out_dir="enrichment"
    ):  # pragma: no cover

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
                pylab.figure()
                ke.scatterplot(direction)
                pylab.tight_layout()
                pylab.savefig(out_dir / f"kegg_{compa}_{direction}.pdf")
                pylab.savefig(out_dir / f"kegg_{compa}_{direction}.png")

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
            try:
                ws.autofilter(0, 0, df.shape[0], df.shape[1] - 1)
            except:
                logger.warning("Fixme")

    def get_gene_lists(self, annot_col="index", Nmax=None):  # pragma: no cover

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

    def _get_design(self, design_file, group, palette):
        """Import design from a table file and add color groups following the
        groups defined in the column 'group' of the table file.
        """
        df = pd.read_csv(design_file, sep=",", index_col=0)
        col_map = dict(zip(df.loc[:, group].unique(), palette))
        df["group_color"] = df.loc[:, group].map(col_map)
        return df

    def _format_plot(self, title="", xlabel="", ylabel="", rotation=0):
        pylab.title(title)
        pylab.xticks(rotation=rotation, ha="right")
        pylab.xlabel(xlabel)
        pylab.ylabel(ylabel)

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

    def plot_count_per_sample(self, fontsize=12, rotation=45):
        """Number of mapped and annotated reads (i.e. counts) per sample. Each color
        for each replicate

        .. plot::
            :include-source:

            from sequana.rnadiff import RNADiffResults
            from sequana import sequana_data

            r = RNADiffResults(sequana_data("rnadiff/rnadiff_onecond_1"))
            r.plot_count_per_sample()

        """
        pylab.clf()
        df = self.counts_raw.sum().rename("total_counts")
        df = pd.concat([self.design_df, df], axis=1)

        pylab.bar(
            df.index,
            df.total_counts / 1000000,
            color=df.group_color,
            lw=1,
            zorder=10,
            ec="k",
            width=0.9,
        )

        pylab.xlabel("Samples", fontsize=fontsize)
        pylab.ylabel("reads (M)", fontsize=fontsize)
        pylab.grid(True, zorder=0)
        pylab.title("Total read count per sample", fontsize=fontsize)
        pylab.xticks(rotation=rotation, ha="right")
        # pylab.xticks(range(N), self.sample_names)
        try:
            pylab.tight_layout()
        except:
            pass

    def plot_percentage_null_read_counts(self):
        """Bars represent the percentage of null counts in each samples.  The dashed
        horizontal line represents the percentage of feature counts being equal
        to zero across all samples

        .. plot::
            :include-source:

            from sequana.rnadiff import RNADiffResults
            from sequana import sequana_data

            r = RNADiffResults(sequana_data("rnadiff/rnadiff_onecond_1"))
            r.plot_percentage_null_read_counts()

        """
        pylab.clf()
        # how many null counts ?
        df = (self.counts_raw == 0).sum() / self.counts_raw.shape[0] * 100
        df = df.rename("percent_null")
        df = pd.concat([self.design_df, df], axis=1)

        pylab.bar(
            df.index, df.percent_null, color=df.group_color, ec="k", lw=1, zorder=10
        )

        all_null = (self.counts_raw == 0).all(axis=1).sum() / self.counts_raw.shape[0]

        pylab.axhline(all_null, ls="--", color="black", alpha=0.5)

        pylab.xticks(rotation=45, ha="right")
        pylab.ylabel("Proportion of null counts (%)")
        pylab.grid(True, zorder=0)
        pylab.tight_layout()

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
            df.index = p.df.columns
            df.columns = ["PC1", "PC2", "PC3"]
            df["size"] = [10] * len(df)  # same size for all points ?

            df = pd.concat([df, self.design_df], axis=1)
            df["label"] = df.index
            df["group_color"] = df["condition"]

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
                text="label",
            )
            return fig
        else:
            variance = p.plot(
                n_components=n_components,
                colors=self.design_df.group_color,
                max_features=max_features,
            )

        return variance

    def plot_mds(self, n_components=2, colors=None, clf=True):
        """IN DEV, not functional"""

        from sequana.viz.mds import MDS

        p = MDS(self.counts_vst)  # [self.sample_names])
        # if colors is None:
        #    colors = {}
        #    for sample in self.sample_names:
        #        colors[sample] = self.colors[self.get_cond_from_sample(sample)]
        p.plot(n_components=n_components, colors=self.design_df.group_color, clf=clf)

    def plot_isomap(self, n_components=2, colors=None):
        """IN DEV, not functional"""

        from sequana.viz.isomap import Isomap

        p = Isomap(self.counts_vst)
        # if colors is None:
        #    colors = {}
        #    for sample in self.sample_names:
        #        colors[sample] = self.colors[self.get_cond_from_sample(sample)]
        p.plot(n_components=n_components, colors=self.design_df.group_color)

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
        df = pd.concat([self.design_df, df], axis=1)

        pylab.clf()
        p = pylab.barh(
            df.index,
            df.most_exp_percent,
            color=df.group_color,
            zorder=10,
            lw=1,
            ec="k",
            height=0.9,
        )

        for idx, rect in enumerate(p):
            pylab.text(
                2,  # * rect.get_height(),
                idx,  # rect.get_x() + rect.get_width() / 2.0,
                df.gene_id.iloc[idx],
                ha="center",
                va="center",
                rotation=0,
                zorder=20,
            )

        self._format_plot(
            # title="Counts monopolized by the most expressed gene",
            # xlabel="Sample",
            xlabel="Percent of total reads",
        )
        pylab.tight_layout()

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
            side_colors=list(self.design_df.group_color.unique()),
        )

        # Convert groups into numbers for Dendrogram category
        group_conv = {
            group: i for i, group in enumerate(self.design_df.condition.unique())
        }
        d.category = self.design_df.condition.map(group_conv).to_dict()
        d.plot()

    def plot_boxplot_rawdata(self, fliersize=2, linewidth=2, rotation=0, **kwargs):
        import seaborn as sbn

        ax = sbn.boxplot(
            data=self.counts_raw.clip(1),
            linewidth=linewidth,
            fliersize=fliersize,
            palette=self.design_df.group_color,
            **kwargs,
        )
        pos, labs = pylab.xticks()
        pylab.xticks(pos, labs, rotation=rotation)
        ax.set_ylabel("Counts (raw) in log10 scale")
        ax.set_yscale("log")
        self._format_plot(ylabel="Raw count distribution")
        pylab.tight_layout()

    def plot_boxplot_normeddata(self, fliersize=2, linewidth=2, rotation=0, **kwargs):
        import seaborn as sbn

        ax = sbn.boxplot(
            data=self.counts_norm.clip(1),
            linewidth=linewidth,
            fliersize=fliersize,
            palette=self.design_df.group_color,
            **kwargs,
        )
        pos, labs = pylab.xticks()
        pylab.xticks(pos, labs, rotation=rotation)
        ax.set(yscale="log")
        self._format_plot(ylabel="Normalised count distribution")
        pylab.tight_layout()

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

    def heatmap(self, comp, log2_fc=1, padj=0.05):
        assert comp in self.comparisons.keys()
        from sequana.viz import heatmap

        h = heatmap.Clustermap(
            self.counts_norm.loc[
                self.comparisons[comp]
                .df.query(
                    "(log2FoldChange<-@log2_fc or log2FoldChange>@log2_fc) and padj<@padj"
                )
                .index
            ]
        ).plot()
