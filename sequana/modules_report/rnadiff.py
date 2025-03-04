# coding: utf-8
#
#  This file is part of Sequana software
#
#  Copyright (c) 2020 - Sequana Development Team
#
#  Distributed under the terms of the 3-clause BSD license.
#  The full license is in the LICENSE file, distributed with this software.
#
#  website: https://github.com/sequana/sequana
#  documentation: http://sequana.readthedocs.io
#
##############################################################################
"""Module to write differential regulation analysis report"""
import os

import colorlog

from sequana.lazy import pandas as pd
from sequana.lazy import pylab
from sequana.modules_report.base_module import SequanaBaseModule
from sequana.rnadiff import RNADiffResults
from sequana.utils.datatables_js import DataTable

logger = colorlog.getLogger(__name__)


__all__ = ["RNAdiffModule"]


class RNAdiffModule(SequanaBaseModule):
    """Write HTML report of variant calling. This class takes a csv file
    generated by sequana_variant_filter.
    """

    def __init__(self, folder, gff, output_filename="summary.html", **kwargs):
        """.. rubric:: constructor"""
        super().__init__()
        self.title = "RNAdiff"
        self.folder = folder
        self.independent_module = True
        self.module_command = "--module rnadiff"

        self.rnadiff = RNADiffResults(folder, gff=gff, **kwargs)

        # nice layout for the report. Use import here to not overload
        # import time
        import seaborn

        seaborn.set()

        # any user option here
        self.kwargs = kwargs

        self._count_section = 1
        self.create_main_report_content()
        self.create_individual_reports()
        self.create_command_section()
        self.create_html(output_filename)

        # Fixme not sure this is required to be imported here
        import matplotlib

        matplotlib.rc_file_defaults()

    def create_command_section(self):
        command = self.kwargs.get("command", "")
        self.sections.append({"name": f"{self._count_section} - Info", "anchor": "command", "content": command})
        self._count_section += 1

    def create_individual_reports(self):
        description = """<p>The differential analysis is based on DESeq2. This
tool aim at fitting one linear model per feature. Given the replicates in
condition one and the replicates in condition two, a p-value is computed to
indicate whether the feature (gene) is differentially expressed. Then, all
p-values are corrected for multiple testing.</p>

<p>It may happen that one sample seems unrelated to the rest. For every feature
and every model, Cook's distance is computed. It reflects how the sample matches
the model. A large value indicates an outlier count and p-values are not computed
for that feature.</p>
"""

        self.sections.append(
            {
                "name": f"{self._count_section}. DGE results",
                "anchor": "filters_option",
                "content": description,
            }
        )
        self._count_section += 1

        counter = 1
        for name, comp in self.rnadiff.comparisons.items():
            self.add_individual_report(comp, name, counter)

    def create_main_report_content(self):
        self.sections = list()

        self.summary()
        self.add_plot_count_per_sample()
        self.add_cluster()
        self.add_normalisation()
        self.add_dispersion()
        if len(self.rnadiff.comparisons) > 1 and len(self.rnadiff.comparisons) < 7:
            self.add_upset_plot()

    def summary(self):
        """Add information of filter."""
        Sdefault = self.rnadiff.summary()
        self.rnadiff.log2_fc = 1
        S1 = self.rnadiff.summary()

        # set options
        options = {
            "scrollX": "true",
            "pageLength": 20,
            "scrollCollapse": "true",
            "dom": "",
            "buttons": [],
        }

        if len(S1) > 20:
            options["dom"] = "Bfrtip"

        def get_local_df(Sdata, lfc=0):
            N = len(Sdata)

            if lfc == 0:
                message = "Number of DGE (any FC)"
            else:
                message = f"Number of DGE log2(|FC|) > {lfc} "

            df = pd.DataFrame(
                {
                    "comparison_link": [1] * N,
                    "comparison": Sdata.index.values,
                    "Description": [message] * N,
                    "Down": Sdata["down"].values,
                    "Up": Sdata["up"].values,
                    "Total": Sdata["all"].values,
                }
            )
            df = df[["comparison", "Description", "Down", "Up", "Total", "comparison_link"]]
            df["comparison_link"] = [f"#{name.replace('.', '')}_stats" for name in Sdata.index]
            return df

        dt = DataTable(get_local_df(Sdefault), "dge_default")
        dt.datatable.set_links_to_column("comparison_link", "comparison", new_page=False)
        dt.datatable.datatable_options = options
        js_all = dt.create_javascript_function()
        html = dt.create_datatable(float_format="%d")

        dt = DataTable(get_local_df(S1, lfc=1), "dge_lfc")
        dt.datatable.set_links_to_column("comparison_link", "comparison", new_page=False)
        dt.datatable.datatable_options = options
        js_all += dt.create_javascript_function()
        html += dt.create_datatable(float_format="%d")

        self.sections.append(
            {
                "name": "Summary",
                "anchor": "filters_option",
                "content": f"""<p>Here below is a summary of the Differential Gene
Expression (DGE) analysis. You can find two entries per comparison. The first
one has no filter except for an adjusted p-value of 0.05. The second shows the
expressed genes with a filter of the log2 fold change of 1 (factor 2 in a normal
scale). Clicking on any of the link will lead you to the section of the comparison.
{js_all} {html} </p>""",
            }
        )

    def add_dispersion(self):
        style = "width:65%"

        def dispersion(filename):
            with pylab.ioff():
                pylab.clf()
                self.rnadiff.plot_dispersion()
                pylab.savefig(filename)
                pylab.close()

        html = """<p>dispersion of the fitted data to the model</p>{}<hr>""".format(
            self.create_embedded_png(dispersion, "filename", style=style)
        )
        self.sections.append({"name": f"{self._count_section}. Dispersion", "anchor": "table", "content": html})
        self._count_section += 1

    def add_cluster(self):
        style = "width:45%"

        def dendogram(filename, count_mode="norm", transform="log"):
            with pylab.ioff():
                pylab.clf()
                import warnings

                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    self.rnadiff.plot_dendogram(count_mode=count_mode, transform_method=transform)
                    try:
                        # Here, tight layout is not comaptible with dendogram
                        # Yet, if we call it and plot again, we get the expected layout.
                        # otherwise, colormap is kind of shifted and not as large as the plot itself.
                        pylab.tight_layout()
                    except Exception:
                        # This does not work right now in v0.16.0 but should be used in future
                        pylab.gcf().set_layout_engine("tight")
                    self.rnadiff.plot_dendogram(count_mode=count_mode, transform_method=transform)

                    pylab.savefig(filename)
                    pylab.close()

        # save image in case of
        dendogram(f"{self.folder}/images/dendogram.png")
        image1 = self.create_embedded_png(dendogram, "filename", style=style)

        if os.path.exists(f"{self.folder}/counts/counts_vst_batch.csv"):

            # save image in case of
            dendogram(f"{self.folder}/images/dendogram_vst_batch.png", count_mode="vst_batch", transform="none")

        # html report
        html_dendogram = f"""<p>The following image shows a hierarchical
clustering of the whole sample set. An euclidean distance is computed between
samples. The dendogram itself is built using the <a
href="https://en.wikipedia.org/wiki/Ward%27s_method"> Ward method </a>. The normed data was log-transformed first and at
max 5,000 most variable features were selected.</p>"""

        if os.path.exists(f"{self.folder}/counts/counts_vst_batch.csv"):
            image2 = self.create_embedded_png(
                dendogram, "filename", count_mode="vst_batch", transform="none", style=style
            )
            html_dendogram += f"""<p>A batch effect was included. Here is the corrected image</p>{image1}{image2}<hr>"""
        else:
            html_dendogram += f"""{image1}<hr>"""

        # =========== PCA
        def pca(filename, transform_method="none", count_mode="vst"):
            pylab.ioff()
            pylab.clf()
            variance = self.rnadiff.plot_pca(
                2,
                fontsize=self.kwargs.get("pca_fontsize", 10),
                transform_method=transform_method,
                count_mode=count_mode,
            )
            pylab.savefig(filename)
            pylab.close()

        image1 = self.create_embedded_png(pca, "filename", style=style)
        html_pca = f"""<p>The experiment variability is also represented by a
principal component analysis as shown here below. The two main components are
represented. We expect the ﬁrst principal component (PC1) to
separate samples from the diﬀerent biological conditions, meaning that the biological variability is
the main source of variance in the data. The data used is the normalised count matrix transformed using a VST method
variance stabilization); the first 500 most variable genes were selected. </p>"""

        from plotly import offline

        html_pca_plotly = "<p>Here is a PCA showing the first 3 components in 3D.</p>"

        import plotly.io as pio

        pio.renderers.default = "iframe"  # Ensure it's set to browser
        pio.renderers["browser"].timeout = 60  # Increa

        if os.path.exists(f"{self.folder}/counts/counts_vst_batch.csv"):
            image2 = self.create_embedded_png(pca, "filename", count_mode="vst_batch", style=style)
            html_pca += f"""<p>Note that a batch effect was included. </p>{image1}{image2}<hr>"""
            fig = self.rnadiff.plot_pca(n_components=3, plotly=True, count_mode="vst_batch")
            html_pca_plotly += offline.plot(fig, output_type="div", include_plotlyjs=False)
        else:
            html_pca += f"""{image1}<hr>"""
            fig = self.rnadiff.plot_pca(n_components=3, plotly=True)
            html_pca_plotly += offline.plot(fig, output_type="div", include_plotlyjs=False)

        self.sections.append(
            {
                "name": f"{self._count_section}. Clusterisation",
                "anchor": "table",
                "content": html_dendogram + html_pca + html_pca_plotly,
            }
        )
        self._count_section += 1

    def add_plot_count_per_sample(self):
        style = "width:65%"

        def plotter(filename):
            pylab.ioff()
            pylab.clf()
            self.rnadiff.plot_count_per_sample(rotation=45)
            pylab.savefig(filename)
            pylab.close()

        html1 = """<p>The following image shows the total number of counted reads
for each sample. We expect counts to be similar within conditions. They may be
different across conditions. Note that variation may happen (e.g., different rRNA contamination
levels, library concentrations, etc).<p>{}<hr>""".format(
            self.create_embedded_png(plotter, "filename", style=style)
        )

        def null_counts(filename):
            pylab.ioff()
            pylab.clf()
            self.rnadiff.plot_percentage_null_read_counts()
            pylab.savefig(filename)
            pylab.close()

        html_null = """<p>The next image shows the percentage of features with no
read count in each sample (taken individually). Features with null read counts
in <b>all</b> samples are not
taken into account in the analysis (black dashed line). Fold-change and p-values
will be set to NA in the final results</p> {}<hr>""".format(
            self.create_embedded_png(null_counts, "filename", style=style)
        )

        def count_density(filename):
            pylab.ioff()
            pylab.clf()
            self.rnadiff.plot_density()
            pylab.savefig(filename)
            pylab.close()

        html_density = """<p>In the following figure, we show the distribution
of read counts for each sample (log10 scale). We expect replicates to behave in
a similar fashion. The mode depends on the biological conditions and organism
considered.</p> {}<hr>""".format(
            self.create_embedded_png(count_density, "filename", style=style)
        )

        def best_count(filename):
            pylab.ioff()
            pylab.clf()
            self.rnadiff.plot_most_expressed_features()
            pylab.savefig(filename)
            pylab.close()

        html_feature = """<p>The following figure shows for each sample (left
y-axis) the gene/feature that captures the highest proportion of the reads
considered. This gene/feature is indicated in the right y-axis. This should not impact
the DESeq2 normalization. We expect consistency across samples within a single
condition</p> {}<hr>""".format(
            self.create_embedded_png(best_count, "filename", style=style)
        )

        self.sections.append(
            {
                "name": f"{self._count_section}. Diagnostic plots",
                "anchor": "table_count",
                "content": html1 + html_null + html_density + html_feature,
            }
        )
        self._count_section += 1

    def add_normalisation(self):
        style = "width:45%"

        def rawcount(filename):
            pylab.ioff()
            pylab.clf()
            self.rnadiff.plot_boxplot_rawdata()
            ax = pylab.gca()
            xticklabels = ax.get_xticklabels()
            ax.set_xticklabels(xticklabels, rotation=45, ha="right")
            try:
                pylab.set_engine_layout("tight")
            except:
                pass
            pylab.savefig(filename)
            pylab.close()

        def normedcount(filename):
            pylab.ioff()
            pylab.clf()
            self.rnadiff.plot_boxplot_normeddata()
            ax = pylab.gca()
            xticklabels = ax.get_xticklabels()
            ax.set_xticklabels(xticklabels, rotation=45, ha="right")
            try:
                pylab.set_engine_layout("tight")
            except:
                pass
            pylab.savefig(filename)
            pylab.close()

        html_boxplot = """<p>A normalization of the data is performed to correct
the systematic technical biases due to different counts across samples. The
normalization is performed with DESeq2. It relies on the hypothess that most
features are not differentially expressed. It computes a scaling factor for each
sample. Normalized read counts are obtained by dividing raw read counts by the
scaling factor associated with the sample they belong to.

Boxplots are often used as a qualitative measure of the quality of the normalization process,
as they show how distributions are globally aﬀected during this process. We expect normalization to
stabilize distributions across samples.
The left figure shows the raw counts while the right figure shows the
normalised counts.
</p>"""
        img1 = self.create_embedded_png(rawcount, "filename", style=style)
        img2 = self.create_embedded_png(normedcount, "filename", style=style)

        html_size_factor = "<p>The size factor is a crucial component used to normalize the raw counts of RNA-seqr data to account for differences in sequencing depth and RNA composition across samples. Each sample is assigned a size factor, which represents the relative amount of sequencing depth (or library size) needed to compare expression levels between samples fairly. This factor adjusts the raw counts by scaling them, so that differences in gene expression are due to biological variation rather than technical artifacts. The size factor is computed by taking the median ratio of each gene's count in a sample to its geometric mean across all samples. Significant departures from 1.0 suggest notable differences in sequencing depth across samples. However, extreme values (either much larger or much smaller than 1.0) might indicate potential issues such as technical variability or inconsistencies in sample preparation. </p>"
        options = {
            "scrollX": "false",
            "pageLength": 1,
            "scrollCollapse": "true",
            "dom": "",
            "buttons": [],
        }
        df = pd.read_csv(f"{self.folder}/code/size_factors.csv", index_col=0)
        datatable = DataTable(df.T, "table_size_factor")
        datatable.datatable.datatable_options = options
        js_all = datatable.create_javascript_function()
        table1 = datatable.create_datatable(float_format="%.3e")

        self.sections.append(
            {
                "name": f"{self._count_section}. Normalisation",
                "anchor": "table",
                "content": html_boxplot + img1 + img2 + html_size_factor + js_all + table1 + "</hr>",
            }
        )
        self._count_section += 1

    def add_upset_plot(self):
        style = "width:65%"

        def upsetplot(filename):
            pylab.ioff()
            pylab.clf()
            self.rnadiff.plot_upset()
            pylab.savefig(filename)
            pylab.close()

        html_upsetplot = """<p> Upset plots are an alternative to venn diagrams, easing the visualisation of DEG lists overlap between comparisons."""
        img = self.create_embedded_png(upsetplot, "filename", style=style)

        self.sections.append(
            {
                "name": f"{self._count_section}. Upset plot",
                "anchor": "table",
                "content": html_upsetplot + img + "</hr>",
            }
        )
        self._count_section += 1

    def add_individual_report(self, comp, name, counter):
        style = "width:45%"

        description = """<p>
When the dispersion estimation and model fitting is done, statistical testing is
performed. The distribution of raw p-values computed by the statistical test
is expected to be a mixture of a uniform distribution on [0, 1] and a peak
around 0 corresponding to the diﬀerentially expressed features. This may not
always be the case. </p>"""

        def plot_pvalue_hist(filename):
            pylab.ioff()
            pylab.clf()
            comp.plot_pvalue_hist()
            pylab.savefig(filename)
            pylab.close()

        def plot_padj_hist(filename):
            pylab.ioff()
            pylab.clf()
            comp.plot_padj_hist()
            pylab.savefig(filename)
            pylab.close()

        img1 = self.create_embedded_png(plot_pvalue_hist, "filename", style=style)
        img2 = self.create_embedded_png(plot_padj_hist, "filename", style=style)

        # FIXME. pvalues adjusted are not relevant so commented for now
        img2 = ""

        self.sections.append(
            {
                "name": f"{self._count_section}.{counter}.a pvalue distribution ({name})",
                "anchor": f"dge_summary",
                "content": description + img1 + img2,
            }
        )

        def plot_volcano(filename):
            pylab.ioff()
            pylab.clf()
            comp.plot_volcano()
            pylab.savefig(filename)
            pylab.close()

        html_volcano = """<p>The volcano plot here below shows the diﬀerentially
expressed features with a adjusted p-value below 0.05 (dashed back line).
The volcano plot represents the log10 of the adjusted P
value as a function of the log2 ratio of diﬀerential expression. </p>"""
        # img3 = self.create_embedded_png(plot_volcano, "filename", style=style)
        img3 = ""

        fig = comp.plot_volcano(
            plotly=True,
            annotations=self.rnadiff.annotation,
            hover_name=self.kwargs.get("hover_name", None),
        )

        from plotly import offline

        plotly = offline.plot(fig, output_type="div", include_plotlyjs=False)
        self.sections.append(
            {
                "name": f"{self._count_section}.{counter}.b volcano plots ({name})",
                "anchor": f"{name}_volcano",
                "content": html_volcano + img3 + "<hr>" + plotly,
            }
        )

        # finally, let us add the tables

        df = comp.df.copy()  # .reset_index()

        # here we need to add the annotation if possible
        try:
            df = pd.concat([df, self.rnadiff.annotation.loc[[str(x) for x in comp.df.index]]], axis=1)
        except Exception as err:
            logger.critical(f"Could not add annotation. {err}. annotation skipped")

        df = df.reset_index()

        fold_change = 2 ** df["log2FoldChange"]
        log10padj = -pylab.log10(df["padj"])
        df.insert(df.columns.get_loc("log2FoldChange") + 1, "FoldChange", fold_change)
        df.insert(df.columns.get_loc("padj") + 1, "log10_padj", log10padj)

        try:
            del df["dispGeneEst"]
        except:
            pass

        for x in ["lfcSE", "stat", "dispersion"]:
            try:
                del df[x]
            except:
                pass
        # set options
        options = {
            "scrollX": "true",
            "pageLength": 10,
            "scrollCollapse": "true",
            "dom": "Bfrtip",
            "buttons": ["copy", "csv"],
        }

        idname = name.replace(".", "")

        datatable = DataTable(df, f"{idname}_table_all")
        datatable.datatable.datatable_options = options
        js_all = datatable.create_javascript_function()
        html_tab_all = datatable.create_datatable(float_format="%.3e")

        df_sign = df.query(f"padj<=0.05 and ({comp.l2fc_name}>1 or {comp.l2fc_name}<-1)")

        datatable = DataTable(df_sign, f"{idname}_table_sign")
        datatable.datatable.datatable_options = options
        js_sign = datatable.create_javascript_function()
        html_tab_sign = datatable.create_datatable(float_format="%.3e")

        content = f"""<p>The following tables give all DGE results. The
first table contains all significant genes (adjusted yyp-value below 0.05 and absolute fold change of at least 0.5). The following tables contains all results
without any filtering. Here is a short explanation for each column:
<ul>
<li> baseMean: read count mean over all samples</li>
<li> FoldChange: fold change in natural base</li>
<li> log2FoldChange: log2 Fold Change estimated by the model. Reflects change between the condition versus the reference condition</li>
<li> stat: Wald statistic for the coefficient (contrast) tested</li>
<li> pvalue: raw p-value from statistical test</li>
<li> padj: adjusted pvalue. Used for cutoff at 0.05 </li>
<li> betaConv: convergence of the coefficients of the model </li>
<li> maxCooks: maximum Cook's distance of the feature </li>
<li> outlier: indicates if the feature is an outlier according to Cook's distance
</li>
</ul>
</p>
<h3>Significative only (and |l2fC| > 1) <a id="{idname}_table_sign"></a></h3>
Here below is a subset of the entire table. It contains all genes below adjusted
p-value of 0.05 and absolute log2 fold change above 1.
{js_sign} {html_tab_sign}"""

        if not self.kwargs.get("split_full_table", False):
            content += f"""<h3>All genes<a id="{idname}_table_all"></a></h3>{js_all} {html_tab_all}"""
        else:
            # create an independent pages (lighter to have N pages rather than everything
            # on the same page. Could be useful for eukaryotes.
            link = f"{self.folder}/{name}.html"
            content += f"""<h3>All genes<a id="{idname}_table_all"></a></h3>"""
            content += f'The full list is available <a target="_blank" href="{name}.html">here</a>'
            r = SequanaBaseModule()
            r.title = "List of DGEs"
            r.folder = self.folder
            r.independent_module = True
            r.sections = list()
            r.sections.append(
                {
                    "name": f"{name} comparison",
                    "anchor": "table",
                    "content": f"""<h3>All genes<a id="{idname}_table_all"></a></h3>{js_all} {html_tab_all}""",
                }
            )
            r.create_html(f"{name}.html")

        self.sections.append(
            {
                "name": f"{self._count_section}.{counter}.c {name} Tables ({name})",
                "anchor": f"{name}_stats",
                "content": content,
            }
        )
