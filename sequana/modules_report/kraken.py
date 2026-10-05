# coding: utf-8
#
#  This file is part of Sequana software
#
#  Copyright (c) 2016 - Sequana Development Team
#
#  File author(s):
#      Thomas Cokelaer <thomas.cokelaer@pasteur.fr>
#      Dimitri Desvillechabrol <dimitri.desvillechabrol@pasteur.fr>,
#          <d.desvillechabrol@gmail.com>
#
#  Distributed under the terms of the 3-clause BSD license.
#  The full license is in the LICENSE file, distributed with this software.
#
#  website: https://github.com/sequana/sequana
#  documentation: http://sequana.readthedocs.io
#
##############################################################################
"""Module to write coverage report"""
import os

import colorlog

from sequana import sequana_data
from sequana.lazy import pandas as pd
from sequana.modules_report.base_module import SequanaBaseModule
from sequana.utils.datatables_js import DataTable

logger = colorlog.getLogger(__name__)


class KrakenModule(SequanaBaseModule):
    """Write HTML report of Kraken results"""

    def __init__(self, input_directory, output_filename=None):
        """
        :param input_directory: the directory of the bwa_bam_to_fastq output
        :param output_filename: if not provided, the HTML is not created.

        """
        super().__init__()
        self.title = "Kraken report"
        self.directory = input_directory
        self.create_report_content()
        if output_filename:
            self.create_html(output_filename)

    def create_report_content(self):
        """Generate the sections list to fill the HTML report."""
        self.sections = list()
        self.add_summary_section()
        self.add_table_results_section()
        self.add_table_results_hierarchy_section()
        self.add_diagnostics_section()
        self.add_blast_section()

    def _get_stats(self):
        return pd.read_csv(self.directory / "kraken.csv")

    def _get_df_hierarchy(self):
        # read the kraken summary results and build a grouped dataframe on ranks
        df = pd.read_csv(self.directory / "kraken.csv")
        df = df.fillna("")
        ranks = [
            "kingdom",
            "phylum",
            "class",
            "order",
            "family",
            "genus",
            "species",
            "name",
        ]
        try:
            df = df.groupby(ranks).sum().sort_index(level=["kingdom"])
        except KeyError:
            return df

        # here we do not want the table to be sorted. it is sorted on first column so
        # let us add back the index
        df = df.reset_index()

        # we also move the taxon/count/percentage at the beginning
        df = df[["taxon", "count", "percentage"] + ranks]
        df = df.reset_index()  # this ensure that order is kept and table not sorted

        # now we fill redundant columns with spaces
        def fix_col(X):
            N = len(X)
            Xp = []
            for i, x in enumerate(X):
                if i == 0:
                    Xp.append(x)
                elif x == X[i - 1]:
                    Xp.append("")
                else:
                    Xp.append(x)
            return Xp

        for rank in ranks:
            df[rank] = fix_col(df[rank])

        return df

    def _get_summary_section(self):
        df = self._get_stats()
        if len(df) == 1 and df.iloc[0]["taxon"] == -1:
            pngimage = sequana_data("no_data.jpg")
            extra = "<p> no reads could be identified with the given the database(s)."
            # Use embedded for bundled resources
            png_html = self.png_to_embedded_png(pngimage)
        else:
            pngimage = self.directory / "kraken.png"
            extra = """<p>The following <b>clickable image</b> is a simplified
version (only genus are shown) of an interactive and more detailled version
based on Krona. Finally, note that the unclassified species in the pie plot
may correspond to species not present in the data base or adapters (if not
removed).</p>"""
            # Use linked PNG (much faster than embedding base64)
            # PNG is in kraken/ subdir, HTML is in parent, so use relative path
            rel_dir = str(self.directory).split(os.sep, 1)[1]  # e.g., "TEST/kraken"
            png_path = rel_dir.split(os.sep)[-1] + "/kraken.png"  # "kraken/kraken.png"
            png_html = self.png_to_linked_png(png_path)

        html = """
    <p>Overview of the Taxonomic content of the filtered reads. </p>
    <p>The taxonomic analysis is performed with Kraken (see database name in
the configuration file. The analysis is performed with a Kmer
approach.
The details about the database itself are available in the <a
href="http://sequana.readthedocs.io">Sequana documentation</a>.
The taxonomic analysis should give a good idea of the content of the FastQ
files but should be used as a sanity check. Indeed, species absent
from the database won't be detected leading to false detection (close species
may be detected instead).
Besides, be aware that closely related species may not be classified precisely.
</p>

    {0}
    <div style="text-align:center"><a href="./{1}/kraken.html"> {2} </a></div>
    <br>
""".format(
            extra,
            str(self.directory).split(os.sep, 1)[1],
            png_html,
        )

        return html

    def _get_table_results(self):
        df = self._get_stats()
        url_ncbi = "https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id={}"
        df["links"] = [url_ncbi.format(taxon) for taxon in df["taxon"]]
        datatable = DataTable(df, "kraken", index=False)

        datatable.datatable.set_links_to_column("links", "taxon")

        datatable.datatable.datatable_options = {
            "scrollX": "300px",
            "pageLength": 30,
            "scrollCollapse": "true",
            "dom": "Bfrtip",
            "paging": "false",
            "order": [[2, "desc"]],
            "buttons": ["copy", "csv"],
        }
        js = datatable.create_javascript_function()
        html_tab = datatable.create_datatable(float_format="%.3g")

        html = f"{html_tab} {js}"
        return html

    def _get_table_grouped_results(self):
        df = self._get_df_hierarchy()
        url_ncbi = "https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id={}"
        df["links"] = [url_ncbi.format(taxon) for taxon in df["taxon"]]
        datatable = DataTable(df, "kraken_summary", index=False)
        datatable.datatable.set_links_to_column("links", "taxon")
        datatable.datatable.datatable_options = {
            "scrollX": "300px",
            "pageLength": 30,
            "scrollCollapse": "true",
            "dom": "frtip",
            "paging": "false",
            # "order": [[ 0, "desc"]],
            "buttons": ["copy", "csv"],
        }
        js = datatable.create_javascript_function()
        html_tab = datatable.create_datatable(float_format="%.3g")
        html = """<p>The following table uses the same data as above. However,
        we have here a hierarchical representation grouping taxons per common lineage.
        <b>DO NOT SORT the table or you will loose the hiearchy</b>"""
        html += f"{html_tab} {js}"

        return html

    def add_summary_section(self):
        html = self._get_summary_section()
        self.sections.append({"name": "Taxonomic content", "anchor": "kraken", "content": html})

    def add_table_results_section(self):
        html = self._get_table_results()
        self.sections.append(
            {
                "name": "Taxonomic Classification",
                "anchor": "classification",
                "content": html,
            }
        )

    def add_table_results_hierarchy_section(self):
        html = self._get_table_grouped_results()
        self.sections.append(
            {
                "name": "Taxonomic Classification (grouped by lineage)",
                "anchor": "grouped_classification",
                "content": html,
            }
        )

    def add_diagnostics_section(self):
        html = """<p>Unclassified read may happen either because the databases used
are not covering the diversity of the sequencing runs, or because reads are of
poor quality. For instance, they may be too short. Here below are information
concerning the read length of unclassified reads. Here below C stands for
classified and U for unclassified reads.</p><div>"""
        # Calculate relative path to kraken subdir (same for all PNGs)
        rel_dir = str(self.directory).split(os.sep, 1)[1]  # e.g., "TEST/kraken"
        subdir = rel_dir.split(os.sep)[-1]  # "kraken"

        pngimage = self.directory / "boxplot_read_length.png"
        if os.path.exists(pngimage):
            html += self.png_to_linked_png(f"{subdir}/boxplot_read_length.png")

        pngimage = self.directory / "hist_read_length.png"
        if os.path.exists(pngimage):
            html += self.png_to_linked_png(f"{subdir}/hist_read_length.png")
        html += "</div>"
        self.sections.append(
            {
                "name": "Diagnostics for unclassified reads",
                "anchor": "diag",
                "content": html,
            }
        )

    def add_blast_section(self):
        logger.debug(f"add_blast_section: directory={self.directory}")
        krona_path = f"{self.directory}/../blast/text.krona.html"
        if os.path.exists(krona_path):
            krona = """ You can visualise the results in the <a href="blast/text.krona.html"> Krona page</a>."""
            logger.debug(f"Krona found: {krona_path}")
        else:
            krona = ""
            logger.debug(f"Krona not found: {krona_path}")

        filename = f"{self.directory}/../blast/blast_summary.csv"
        logger.debug(f"Looking for blast CSV: {filename}, exists={os.path.exists(filename)}")
        try:
            if os.path.exists(filename):
                df = pd.read_csv(filename, sep=",")
                logger.debug(f"Loaded blast CSV, shape={df.shape}")
                datatable = DataTable(df, "blast_summary", index=False)
                datatable.datatable.datatable_options = {
                    "scrollX": "300px",
                    "pageLength": 30,
                    "scrollCollapse": "true",
                    "dom": "frtip",
                    "paging": "false",
                    "buttons": ["copy", "csv"],
                }
                js = datatable.create_javascript_function()
                html_tab = datatable.create_datatable(float_format="%.3g")
                html = f"""<p>Here below you can find the summary of unclassified reads blasted on a local BLAST database. {krona}</p>
{html_tab}
{js}
<hr>"""
                logger.debug(f"Created blast table, html length={len(html)}")
            else:
                html = f"{krona}"
                logger.debug(f"Blast CSV not found, using krona only")
        except Exception as e:
            logger.warning(f"Could not create blast table: {e}")
            logger.exception("Full traceback:")
            html = f"{krona}"

        self.sections.append(
            {
                "name": "Blast of unclassified reads overview",
                "anchor": "blast",
                "content": html,
            }
        )
