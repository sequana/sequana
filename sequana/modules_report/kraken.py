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
import glob
import io

from sequana.modules_report.base_module import SequanaBaseModule
from sequana.utils import config

from sequana.lazy import pandas as pd
from sequana.lazy import pylab
from sequana import logger, sequana_data

from sequana.utils.datatables_js import DataTable


class KrakenModule(SequanaBaseModule):
    """ Write HTML report of Kraken results"""
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
        """ Generate the sections list to fill the HTML report.
        """
        self.sections = list()
        self.add_summary_section()

    def _get_stats(self):
        df = pd.read_csv(self.directory + os.sep + "kraken.csv")
        return df

    def _get_summary_section(self):


        df = self._get_stats()
        if len(df) == 1 and df.iloc[0]['taxon'] == -1:
            pngimage = sequana_data("no_data.jpg")
            extra = "<p> no reads could be identified with the given the database(s)."
        else:
            pngimage = self.directory + os.sep + "kraken.png"
            extra = """<p>The following <b>clickable image</b> is a simplified 
version (only genus are shown) of an interactive and more detailled version 
based on Krona. Finally, note that the unclassified species in the pie plot 
may correspond to species not present in the data base or adapters (if not 
removed).</p>"""

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
""".format(extra, self.directory.split(os.sep, 1)[1], 
            self.png_to_embedded_png(pngimage))

        datatable = DataTable(df, "kraken", index=False)
        # add links
        if "ena" in df.columns:
            urlena = "http://www.ebi.ac.uk/ena/data/view/"
            datatable.datatable.set_links_to_column("ena",
                [urlena + this for this in df['ena']])
        datatable.datatable.datatable_options = {
            'scrollX': '300px',
            'pageLength': 30,
            'scrollCollapse': 'true',
            'dom': 'irtpB',
            "paging": "false",
            "order": [[ 2, "desc"]],
            'buttons': ['copy', 'csv']}
        js = datatable.create_javascript_function()
        html_tab = datatable.create_datatable(float_format='%.3g')

        html += "{} {}".format(html_tab, js)
        """# Rounding and convert in string to avoid exp notation
        df['percentage']  = df['percentage'].apply(lambda x: str(round(x,4)))
        #self.jinja['kraken_json'] = df.to_json()"""

        return html

    def add_summary_section(self):
        html = self._get_summary_section()
        self.sections.append({
          "name": "Taxonomic content",
          "anchor": "kraken",
          "content": html
        })

