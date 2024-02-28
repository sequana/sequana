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
"""Module to write enrichment report"""
import os
import sys
from pathlib import Path

import colorlog

from sequana.lazy import pandas as pd
from sequana.lazy import pylab
from sequana.modules_report.base_module import SequanaBaseModule
from sequana.utils import config
from sequana.utils.datatables_js import DataTable

logger = colorlog.getLogger(__name__)


class ModuleUniprotEnrichment(SequanaBaseModule):
    """Write HTML report of variant calling. This class takes a csv file
    generated by sequana_variant_filter.
    """

    def __init__(
        self,
        gene_lists,
        summary,
        enrichment_params={
            "padj": 0.05,
            "log2_fc": 3,
            "max_entries": 3000,
            "nmax": 50,
            "plot_logx": True,
        },
        command="",
        ontologies=["MF", "BP", "CC"],
    ):
        """.. rubric:: constructor"""
        super().__init__()
        self.title = "UniProt Enrichment (GO terms)"
        self.taxon_id = summary["taxon"]
        self.taxon_name = summary["taxon_name"]
        self.compa = summary["name"]

        self.command = command
        self.gene_lists = gene_lists
        self.enrichment_params = enrichment_params
        self.nmax = enrichment_params.get("nmax", 50)

        # compute the enrichment here once for all, This may take time
        from sequana import logger

        logger.setLevel("INFO")

        logger.info(" === Module UniProtEnrichment. === ")
        self.ontologies = ontologies

        self.create_report_content()
        print(f"{config.output_dir}/{self.compa}_enrichment.html")
        # config.output_dir used internally
        self.create_html(f"{self.compa}_enrichment.html")

    def create_report_content(self):
        self.sections = list()
        self.summary()
        self.add_go()
        self.sections.append({"name": "5 - Info", "anchor": "command", "content": self.command})

    def summary(self):
        """Add information."""

        total_up = len(self.gene_lists["up"])
        total_down = len(self.gene_lists["down"])
        total = total_up + total_down
        log2fc = self.enrichment_params["log2_fc"]

        # create html table for taxon information
        _taxon_id = self.taxon_id
        _taxon_name = self.taxon_name

        js = ""
        html_table = ""
        self.sections.append(
            {
                "name": "1 - Summary",
                "anchor": "filters_option",
                "content": f"""

<p>In the following sections, you will find the GO
terms enrichment. The input data for those analysis is the output of the RNADiff
analysis where adjusted p-values above 0.05 are excluded. Moreover, we removed
candidates with log2 fold change below {log2fc}. Using these filters, the list of
differentially expressed genes is made of {total_up} up and {total_down} down genes (total {total})</p>
<p> In the following plots you can find the first GO terms that are enriched, keeping a
maximum of {self.nmax} identifiers. </p>

<p>The taxon used is {_taxon_name} (ID {_taxon_id}).<br>

""",
            }
        )

    def add_go(self):
        # somehow, logger used here and in add_kegg must be global. If you call
        # add_go and then add_kegg, logger becomes an unbound local variable.
        # https://stackoverflow.com/questions/10851906/python-3-unboundlocalerror-local-variable-referenced-before-assignment
        # global logger
        level = logger.level
        logger.setLevel(level)

        html_intro = """
<p>Here below is a set of plots showing the enriched GO terms using the down
regulated genes only, and then the up-regulated genes only. When possible a
graph of the found GO terms is provided. MF stands for molecular
function, CC for cellular components and BP for biological process.</p>
        </div>
        """

        html = self._get_enrichment("down")
        self.sections.append(
            {
                "name": "2 - Enriched GO terms (Down cases)",
                "anchor": "go_down",
                "content": html,
            }
        )

        html = self._get_enrichment("up")
        self.sections.append(
            {
                "name": "3 - Enriched GO terms (Up cases)",
                "anchor": "go_up",
                "content": html,
            }
        )

        html = self._get_enrichment("all")
        self.sections.append(
            {
                "name": "4 - Enriched GO terms (All cases)",
                "anchor": "go_all",
                "content": html,
            }
        )

    # a utility function to create the proper html table
    def get_html_table(self, this_df, identifier):
        df = this_df.copy()

        links = []
        for x in df["Term"]:
            links.append(f"https://www.ebi.ac.uk/QuickGO/term/{x}")
        df["links"] = links

        # remove non-informative or redundant fields
        df = df.drop(
            ["Gene_set", "Unnamed: 0", "id", "abs_log2_fold_enrichment"],
            errors="ignore",
            axis=1,
        )

        # we move the "Genes" to the end
        df = df[[x for x in df.columns if x != "Genes"] + ["Genes"]]

        first_col = df.pop("Term")
        df.insert(0, "Term", first_col)
        df = df.sort_values(by="fold_enrichment", ascending=False)

        datatable = DataTable(pd.DataFrame(df), identifier)
        datatable.datatable.set_links_to_column("links", "Term")
        datatable.datatable.datatable_options = {
            "scrollX": "true",
            "pageLength": 10,
            "scrollCollapse": "true",
            "dom": "Bfrtip",
            "buttons": ["copy", "csv"],
        }
        js = datatable.create_javascript_function()
        html_table = datatable.create_datatable(float_format="%E")
        return js + html_table

    def _get_enrichment(self, category):
        # category is in down/up/all

        style = "width:95%"

        html = ""

        filenames = []

        for ontology in self.ontologies:
            filename = f"{config.output_dir}/plot_{category}_{ontology}.png"

            if os.path.exists(filename):
                df = pd.read_csv(f"{config.output_dir}/{category}_{ontology}.csv")
                image = self.png_to_embedded_png(
                    filename,
                    style=style,
                )

                n_go_terms = len(df)
                html += f"""
<h3>{category.title()} - {ontology}</h3>
<p>For {ontology}, we found {n_go_terms} go terms.
Showing {self.nmax} here below (at most). The full list is downlodable from the CSV
 file hereafter.</p> {image} <br>"""
                html += self.get_html_table(df, f"GO_table_{category}_{ontology}")
                filenames.append(Path(f"chart_{category}_{ontology}.png"))
            else:
                html += f"""
<h4>{category.title()} - {ontology}</h4><p>For {ontology} case, we found 0
enriched go terms. </p><br>"""

        foto = self.add_fotorama(filenames, width=1000)
        html += f"<h4>Charts {category} -- </h4> {foto}"

        return html
