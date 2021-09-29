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
import subprocess
from pathlib import Path
from itertools import combinations

from jinja2 import Environment, PackageLoader
import seaborn as sns

from sequana.lazy import pandas as pd
from sequana.lazy import pylab
from sequana.lazy import numpy as np
from sequana.gff3 import GFF3
from sequana.viz import Volcano
from sequana.enrichment import PantherEnrichment
from sequana.enrichment import KEGGPathwayEnrichment
from sequana.featurecounts import FeatureCount


import colorlog

logger = colorlog.getLogger(__name__)


__all__ = ["RNADiffAnalysis", "RNADiffResults", "RNADiffTable", "RNADesign"]


def strip(text):
    try:
        return text.strip()
    except AttributeError:
        return text


class RNADesign:
    """Simple RNA design handler"""

    def __init__(self, filename, sep=r"\s*,\s*", reference=None):
        self.filename = filename
        # \s to strip the white spaces
        self.df = pd.read_csv(
            filename, sep=sep, engine="python", comment="#", dtype={"label": str}
        )
        if reference and reference not in self.conditions:
            raise ValueError(
                f"{reference} condition (for the reference) not found in the conditions of your design file"
            )
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

    def keep_conditions(self, conditions):
        self.df = self.df.query("condition in @conditions")


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
    :param keep_all_conditions: if user set comparisons, it means will only want
        to include some comparisons and therefore their conditions. Yet,
        sometimes, you may still want to keep all conditions in the diffential
        analysis. If some set this flag to True.
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
    :param sep_counts: The separator used in the input count file.
    :param sep_design: The separator used in the input design file.

    This class reads a :class:`sequana.featurecounts.`

    ::

        r = rnadiff.RNADiffAnalysis("counts.csv", "design.csv",
                condition="condition", comparisons=[(("A", "B"), ('A', "C")],


    For developers: the rnadiff_template.R script behind the scene expects those
    attributes to be found in the RNADiffAnalysis class: counts_filename,
    design_filename, fit_type, fonction, comparison_str, independent_filtering,
    cooks_cutoff, code_dir, outdir, counts_dir, beta_prior, threads

    """

    _template_file = "rnadiff_light_template.R"
    _template_env = Environment(loader=PackageLoader("sequana", "resources/templates"))
    template = _template_env.get_template(_template_file)

    def __init__(
        self,
        counts_file,
        design_file,
        condition,
        keep_all_conditions=False,
        reference=None,
        comparisons=None,
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
        sep_design=r"\s*,\s*",
        minimum_mean_reads_per_gene=0,
    ):

        # if set, we can filter genes that have low counts (on average)
        self.minimum_mean_reads_per_gene = minimum_mean_reads_per_gene

        # define some output directory and create them
        self.outdir = Path(outdir)
        self.counts_dir = self.outdir / "counts"
        self.code_dir = self.outdir / "code"

        self.outdir.mkdir(exist_ok=True)
        self.code_dir.mkdir(exist_ok=True)
        self.counts_dir.mkdir(exist_ok=True)

        self.usr_counts = counts_file

        self.counts_filename = self.code_dir / "counts.csv"

        # Read and check the design file. Filtering if comparisons is provided
        self.design = RNADesign(design_file, sep=sep_design, reference=reference)
        self.comparisons = comparisons if comparisons else self.design.comparisons


        _conditions = {x for comp in self.comparisons for x in comp}
        if not keep_all_conditions:
            self.design.keep_conditions(_conditions)
        logger.info(f"Conditions that are going to be included: ")
        for x in self.design.conditions:
            logger.info(f" - {x}")
        # we do not sort the design but the user order. Important for plotting
        self.design = self.design.df.set_index("label")

        # save the design file keeping track of its name
        self.design_filename = self.code_dir / "design.csv"
        self.design.to_csv(self.design_filename)

        # Reads and check the count file
        self.counts = self.check_and_save_input_tables(sep_counts)

        # the name of the condition in the design file
        self.condition = condition
        self.check_condition()

        # check comparisons and print information
        self.check_comparisons()

        logger.info(f"Comparisons to be included:")
        for x in self.comparisons:
            logger.info(f" - {x}")
        self.comparisons_str = (
            f"list({', '.join(['c' + str(x) for x in self.comparisons])})"
        )

        # For DeSeq2
        self.batch = batch
        self.model = f"~{batch + '+' + condition if batch else condition}"
        self.fit_type = fit_type
        self.beta_prior = "TRUE" if beta_prior else "FALSE"
        self.independent_filtering = "TRUE" if independent_filtering else "FALSE"
        self.cooks_cutoff = cooks_cutoff if cooks_cutoff else "TRUE"

        # for metadata
        self.gff = gff
        self.fc_feature = fc_feature
        self.fc_attribute = fc_attribute
        self.annot_cols = annot_cols
        self.threads = threads

        # sanity check for the R scripts:
        for attr in (
            "counts_filename",
            "design_filename",
            "fit_type",
            "comparisons_str",
            "independent_filtering",
            "cooks_cutoff",
            "code_dir",
            "outdir",
            "counts_dir",
            "beta_prior",
            "threads",
        ):
            try:
                getattr(self, attr)
            except AttributeError as err:
                logger.error(
                    f"Attribute {attr} missing in the RNADiffAnalysis class. cannot go further"
                )
                raise Exception(err)

    def __repr__(self):
        info = f"RNADiffAnalysis object:\n\
- {self.counts.shape[1]} samples.\n\
- {len(self.comparisons)} comparisons.\n\n\
Counts overview:\n\
{self.counts.head()}\n\n\
Design overview:\n\
{self.design.head()}"

        return info

    def check_and_save_input_tables(self, sep_counts):

        # input may be an existing rnadiff.csv file create with FeatureCount
        # class, or (if it fails) a feature count file (tabulated with
        # Chr/Start/Geneid columns)
        try:
            # low memory set to False to avoid warnings (see pandas.read_csv doc)
            counts = pd.read_csv(
                self.usr_counts,
                sep=sep_counts,
                index_col="Geneid",
                comment="#",
                low_memory=False,
            )
        except ValueError:
            counts = FeatureCount(self.usr_counts).df

        if self.minimum_mean_reads_per_gene > 0:
            logger.info(f"{len(counts)} annotated feature to be processed")
        counts = counts[counts.mean(axis=1) > self.minimum_mean_reads_per_gene]
        if self.minimum_mean_reads_per_gene > 0:
            logger.info(f"Keeping {len(counts)} features after removing low "
                f"counts below {self.minimum_mean_reads_per_gene} on average")

        # filter count based on the design and the comparisons provide in the
        # constructor so that columns match as expected by a DESeq2 analysis.

        counts = counts[self.design.index]

        # Save this sub count file

        counts.to_csv(self.counts_filename)

        return counts

    def check_condition(self):
        # set condition of the statistical model

        columns = ",".join(self.design.columns)
        if self.condition not in columns:
            logger.error(
                f"""Your condition named '{condition}' is expected to
be in the header of your design file but was not found. Candidates are:
    {columns}"""
            )
            sys.exit(1)

    def check_comparisons(self):
        # let us check the consistenct of the design and comparisons
        valid_conditions = ",".join(set(self.design[self.condition].values))
        for item in [x for y in self.comparisons for x in y]:
            if item not in self.design[self.condition].values:
                logger.error(
                    f"""{item} not found in the design. Fix the design
or comparisons. possible values are {valid_conditions}"""
                )
                sys.exit(1)

    def run(self):
        """Create outdir and a DESeq2 script from template for analysis. Then execute
        this script.

        :return: a :class:`RNADiffResults` instance
        """
        logger.info("Running DESeq2 analysis. Please wait")

        rnadiff_script = self.code_dir / "rnadiff_light.R"

        with open(rnadiff_script, "w") as f:
            f.write(RNADiffAnalysis.template.render(self.__dict__))

        logger.info("Starting differential analysis with DESeq2...")

        # capture_output is valid for py3.7 and above so we will use
        # stdout/stderr to be back compatible with py3.6
        # Capture rnadiff output
        p = subprocess.Popen(
            f"Rscript {rnadiff_script}",
            shell=True,
            universal_newlines=True,
            stderr=subprocess.PIPE,
            stdout=subprocess.PIPE,
        )
        p.wait()
        stdout, stderr = p.stdout.read(), p.stderr.read()

        # Capture rnadiff output, Unfortunately, R code mixes stdout/stderr
        # FIXME
        with open(self.code_dir / "rnadiff.err", "w") as f:
            f.write(stderr)
        with open(self.code_dir / "rnadiff.out", "w") as f:
            f.write(stdout)

        logger.info("DGE analysis done.")

        results = RNADiffResults(
            self.outdir,
            condition=self.condition,
            gff=self.gff,
            fc_feature=self.fc_feature,
            fc_attribute=self.fc_attribute,
            annot_cols=self.annot_cols,
        )
        return results


class RNADiffTable:
    def __init__(self, path, alpha=0.05, log2_fc=0, sep=",", condition="condition"):
        """A representation of the results of a single rnadiff comparison

        Expect to find output of RNADiffAnalysis file named after condt1_vs_cond2_degs_DESeq2.csv

        ::

            from sequana.rnadiff import RNADiffTable
            RNADiffTable("A_vs_B_degs_DESeq2.csv")


        """
        self.path = Path(path)
        self.name = self.path.stem.replace("_degs_DESeq2", "").replace("-", "_")

        self._alpha = alpha
        self._log2_fc = log2_fc

        self.df = pd.read_csv(self.path, index_col=0, sep=sep)
        self.df.padj[self.df.padj == 0] = 1e-50
        self.condition = condition

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

        fc_filt = self.df["log2FoldChange"].abs() < self._log2_fc
        fdr_filt = self.df["padj"] > self._alpha
        outliers = self.df["padj"].isna()

        filt_df = self.df.copy()
        filt_df[fc_filt.values | fdr_filt.values | outliers] = np.NaN
        return filt_df

    def set_gene_lists(self):

        only_drgs_df = self.filt_df.dropna(how="all")

        self.gene_lists = {
            "up": list(only_drgs_df.query("log2FoldChange > 0").index),
            "down": list(only_drgs_df.query("log2FoldChange < 0").index),
            "all": list(only_drgs_df.index),
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
        annotations=None,
        hover_name=None,
    ):
        """

        .. plot::
            :include-source:

            from sequana.rnadiff import RNADiffResults
            from sequana import sequana_data

            r = RNADiffResults(sequana_data("rnadiff/", "doc"))
            r.comparisons["A_vs_B"].plot_volcano()

        """

        if plotly:
            from plotly import express as px

            df = self.df.copy()

            if annotations is not None:
                try:
                    df = pd.concat([df, annotations.annotation], axis=1)
                except Exception as err:
                    logger.warning(
                        f"Could not merge rnadiff table with annotation. Full error is: {err}"
                    )
            df["log_adj_pvalue"] = -pylab.log10(df.padj)
            df["significance"] = [
                "<{}".format(padj) if x else ">={}".format(padj) for x in df.padj < padj
            ]

            if hover_name is not None:
                if hover_name not in df.columns:
                    logger.warning(
                        f"hover_name {hover_name} not in the GFF attributes. Switching to automatic choice"
                    )
                    hover_name = None
            if hover_name is None:
                for name in ["Name", "gene_name", "gene_id", "locus_tag", "ID"]:
                    if name in df.columns:
                        hover_name = name

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

        d1 = self.df.query("padj>@padj")
        d2 = self.df.query("padj<=@padj")
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
            df[self.condition] = [
                self.get_cond_from_sample(sample) for sample in self.sample_names
            ]
            fig = px.scatter_3d(
                df,
                x="PC1",
                y="PC2",
                z="PC3",
                color=self.condition,
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
        gff=None,
        fc_attribute=None,
        fc_feature=None,
        pattern="*vs*_degs_DESeq2.csv",
        alpha=0.05,
        log2_fc=0,
        palette=sns.color_palette(desat=0.6),
        condition="condition",
        annot_cols=None,
        # annot_cols=["ID", "Name", "gene_biotype"],
        **kwargs,
    ):
        """

        :rnadiff_folder: a valid rnadiff folder created by :class:`RNADiffAnalysis`

        ::

            RNADiffResults("rnadiff/")


        """
        self.path = Path(rnadiff_folder)
        self.files = [x for x in self.path.glob(pattern)]

        self.counts_raw = pd.read_csv(
            self.path / "counts" / "counts_raw.csv", index_col=0, sep=","
        )
        self.counts_raw.sort_index(axis=1, inplace=True)

        self.counts_norm = pd.read_csv(
            self.path / "counts" / "counts_normed.csv", index_col=0, sep=","
        )
        self.counts_norm.sort_index(axis=1, inplace=True)

        self.counts_vst = pd.read_csv(
            self.path / "counts" / "counts_vst_norm.csv", index_col=0, sep=","
        )
        self.counts_vst.sort_index(axis=1, inplace=True)

        self.dds_stats = pd.read_csv(
            self.path / "code" / "overall_dds.csv", index_col=0, sep=","
        )
        self.condition = condition

        design_file = f"{rnadiff_folder}/code/design.csv"
        self.design_df = self._get_design(
            design_file, condition=self.condition, palette=palette
        )

        # optional annotation
        self.fc_attribute = fc_attribute
        self.fc_feature = fc_feature
        self.annot_cols = annot_cols
        if gff:
            if fc_feature is None or fc_attribute is None:
                logger.warning(
                    "Since you provided a GFF file you must provide the feature and attribute to be used."
                )
            self.annotation = self.read_annot(gff)
        else:
            self.annotation = pd.read_csv(
                self.path / "rnadiff.csv", index_col=0, header=[0, 1]
            )
            self.annotation = self.annotation["annotation"]

        # some filtering attributes
        self._alpha = alpha
        self._log2_fc = log2_fc

        self.comparisons = self.import_tables()

        self.df = self._get_total_df()
        self.filt_df = self._get_total_df(filtered=True)

        self.fontsize = kwargs.get("fontsize", 12)
        self.xticks_fontsize = kwargs.get("xticks_fontsize", 12)

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
        logger.warning("DEPRECATED DO NOT USE read_csv from RNADiffResults")
        self.df = pd.read_csv(filename, index_col=0, header=[0, 1])

    def import_tables(self):
        from easydev import AttrDict

        data = {
            compa.stem.replace("_degs_DESeq2", "").replace("-", "_"): RNADiffTable(
                compa,
                alpha=self._alpha,
                log2_fc=self._log2_fc,
                condition=self.condition,
                # gff=self.annotation.annotation,
            )
            for compa in self.files
        }

        return AttrDict(**data)

    def read_annot(self, gff):
        """Get a properly formatted dataframe from the gff.

        :param gff: a input GFF filename or an existing instance of GFF3
        """

        # if gff is already instanciate, we can just make a copy otherwise
        # we read it indeed.
        if not hasattr(gff, "df"):
            gff = GFF3(gff)

        if self.annot_cols is None:
            lol = [
                list(x.keys())
                for x in gff.df.query("type==@self.fc_feature")["attributes"].values
            ]
            annot_cols = sorted(list(set([x for item in lol for x in item])))
        else:
            annot_cols = self.annot_cols

        df = gff.df.query("type == @self.fc_feature").loc[:, annot_cols]
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

    def get_gene_lists(
        self, annot_col="index", Nmax=None, dropna=False
    ):  # pragma: no cover

        gene_lists_dict = {}

        for compa in self.comparisons.keys():
            df = self.df.loc[:, [compa]].copy()
            df = df.droplevel(0, axis=1)

            # Let us add the annotation columns
            df = pd.concat([df, self.annotation.loc[df.index]], axis=1)

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

            if dropna:
                for direction in gene_lists_dict[compa]:
                    gl = gene_lists_dict[compa][direction]
                    if not gl:
                        continue
                    perc_unannotated = gl.count(None) / len(gl) * 100
                    logger.warning(
                        f"{compa} {direction}: Removing {perc_unannotated:.0f}% of the data for enrichment analysis due to missing identifiers in annotation."
                    )
                    gene_lists_dict[compa][direction] = [x for x in gl if x != None]

        return gene_lists_dict

    def _get_design(self, design_file, condition, palette):
        """Import design from a table file and add color groups following the
        groups defined in the column 'condition' of the table file.
        """
        design = RNADesign(design_file)
        df = design.df.set_index("label")

        if len(design.conditions) > len(palette):
            palette = sns.color_palette("deep", n_colors=len(design.conditions))

        col_map = dict(zip(df.loc[:, condition].unique(), palette))

        df["group_color"] = df.loc[:, condition].map(col_map)
        return df

    def _format_plot(self, title="", xlabel="", ylabel="", rotation=0, fontsize=None):
        pylab.title(title)
        pylab.xticks(rotation=rotation, ha="right", fontsize=fontsize)
        pylab.xlabel(xlabel)
        pylab.ylabel(ylabel)

    def get_specific_commons(self, direction, compas=None, annot_col="index"):
        """Extract gene lists for all comparisons.

        Genes are common (but specific, ie a gene appear only appears in the
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

    def _set_figsize(self, height=5, width=8):
        pylab.figure()
        fig = pylab.gcf()
        fig.set_figheight(height)
        fig.set_figwidth(width)

    def plot_count_per_sample(self, fontsize=None, rotation=45, xticks_fontsize=None):
        """Number of mapped and annotated reads (i.e. counts) per sample. Each color
        for each replicate

        .. plot::
            :include-source:

            from sequana.rnadiff import RNADiffResults
            from sequana import sequana_data

            r = RNADiffResults(sequana_data("rnadiff/", "doc"))
            r.plot_count_per_sample()

        """
        self._set_figsize()
        if fontsize is None:
            fontsize = self.fontsize
        if xticks_fontsize is None:
            xticks_fontsize = self.xticks_fontsize

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
        pylab.xticks(rotation=rotation, ha="right", fontsize=xticks_fontsize)
        try:
            pylab.tight_layout()
        except:
            pass

    def plot_percentage_null_read_counts(self, fontsize=None, xticks_fontsize=None):
        """Bars represent the percentage of null counts in each samples.  The dashed
        horizontal line represents the percentage of feature counts being equal
        to zero across all samples

        .. plot::
            :include-source:

            from sequana.rnadiff import RNADiffResults
            from sequana import sequana_data

            r = RNADiffResults(sequana_data("rnadiff/", "doc"))
            r.plot_percentage_null_read_counts()

        """
        self._set_figsize()
        if fontsize is None:
            fontsize = self.fontsize
        if xticks_fontsize is None:
            xticks_fontsize = self.xticks_fontsize

        # how many null counts ?
        df = (self.counts_raw == 0).sum() / self.counts_raw.shape[0] * 100
        df = df.rename("percent_null")
        df = pd.concat([self.design_df, df], axis=1)

        pylab.bar(
            df.index, df.percent_null, color=df.group_color, ec="k", lw=1, zorder=10
        )

        all_null = (self.counts_raw == 0).all(axis=1).sum() / self.counts_raw.shape[0]

        pylab.axhline(all_null, ls="--", color="black", alpha=0.5)

        pylab.xticks(rotation=45, ha="right", fontsize=xticks_fontsize)
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
        fontsize=10,
        adjust=True
    ):

        """

        .. plot::
            :include-source:

            from sequana.rnadiff import RNADiffResults
            from sequana import sequana_data

            r = RNADiffResults(sequana_data("rnadiff/", "doc"))

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
            df["group_color"] = df[self.condition]

            # plotly uses 10 colors by default. Here we cope with the case
            # of having more than 10 conditions
            colors = None
            try:
                if len(set(self.design_df.condition.values)):
                    colors = sns.color_palette("deep", n_colors=13)

                    colors = px.colors.qualitative.Light24
            except Exception as err:
                logger.warning("Could not determine number of conditions")

            fig = px.scatter_3d(
                df,
                x="PC1",
                y="PC2",
                z="PC3",
                color="group_color",
                color_discrete_sequence=colors,
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
                fontsize=fontsize,
                adjust=adjust
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

    def plot_feature_most_present(self, fontsize=None, xticks_fontsize=None):
        """"""
        if fontsize is None:
            fontsize = self.fontsize
        if xticks_fontsize is None:
            xticks_fontsize = self.xticks_fontsize

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
        pylab.yticks(fontsize=xticks_fontsize)

        for idx, rect in enumerate(p):
            pylab.text(
                2,  # * rect.get_height(),
                idx,  # rect.get_x() + rect.get_width() / 2.0,
                df.gene_id.iloc[idx],
                ha="center",
                va="center",
                rotation=0,
                zorder=20,
                fontsize=xticks_fontsize
            )

        self._format_plot(
            # title="Counts monopolized by the most expressed gene",
            # xlabel="Sample",
            xlabel="Percent of total reads",
            fontsize=xticks_fontsize
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
            group: i for i, group in enumerate(self.design_df[self.condition].unique())
        }
        d.category = self.design_df[self.condition].map(group_conv).to_dict()
        d.plot()

    def plot_boxplot_rawdata(self, fliersize=2, linewidth=2, rotation=0,
            fontsize=None, xticks_fontsize=None, **kwargs):

        import seaborn as sbn
        if fontsize is None:
            fontsize = self.fontsize
        if xticks_fontsize is None:
            xticks_fontsize = self.xticks_fontsize

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
        self._format_plot(ylabel="Raw count distribution", fontsize=xticks_fontsize)
        pylab.tight_layout()

    def plot_boxplot_normeddata(self, fliersize=2, linewidth=2, rotation=0,
        fontsize=None, xticks_fontsize=None, **kwargs):

        import seaborn as sbn
        if fontsize is None:
            fontsize = self.fontsize
        if xticks_fontsize is None:
            xticks_fontsize = self.xticks_fontsize

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

    def _replace_index_with_annotation(self, df, annot):
        # ID is unique but annotation_column may not be complete with NA
        # Let us first get the annotion with index as the data index
        # and one column (the annotation itself)
        dd = self.annotation.loc[df.index][annot]

        # Let us replace the possible NA with the ID
        dd = dd.fillna(dict({(x, x) for x in dd.index}))

        # Now we replace the data index with this annoation
        df.index = dd.values

        return df

    def heatmap_vst_centered_data(
        self,
        comp,
        log2_fc=1,
        padj=0.05,
        xlabel_size=8,
        ylabel_size=12,
        figsize=(10, 15),
        annotation_column=None,
    ):

        assert comp in self.comparisons.keys()
        from sequana.viz import heatmap

        # Select counts based on the log2 fold change and padjusted
        data = self.comparisons[comp].df.query(
            "(log2FoldChange<-@log2_fc or log2FoldChange>@log2_fc) and padj<@padj"
        )
        counts = self.counts_vst.loc[data.index].copy()

        logger.info(f"Using {len(data)} DGE genes")

        # replace the indices with the proper annotation if required.
        if annotation_column:
            data = self._replace_index_with_annotation(data, annotation_column)
            counts.index = data.index

        # finally the plots
        h = heatmap.Clustermap(counts, figsize=figsize, z_score=0, center=0)

        ax = h.plot()
        ax.ax_heatmap.tick_params(labelsize=xlabel_size, axis="x")
        ax.ax_heatmap.tick_params(labelsize=ylabel_size, axis="y")

        return ax
