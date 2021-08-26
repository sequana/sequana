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
import io
import glob

from sequana.modules_report.base_module import SequanaBaseModule

from sequana.lazy import pandas as pd
from sequana.lazy import pylab

import colorlog

logger = colorlog.getLogger(__name__)


from sequana.utils.datatables_js import DataTable


class FastQCModule(SequanaBaseModule):
    """Write HTML report for fastqc.

    Searches for _fastqc.html files

    """

    def __init__(self, output_filename="fastqc.html", pattern="*/*_fastqc.html"):
        """

        :param input:
        :param pattern: we use a glob to search for the relevant files
        """
        super().__init__()
        self.title = "FastQC"
        self.pattern = pattern
        self.create_report_content()
        self.create_html(output_filename)

    def create_report_content(self):
        self.sections = list()
        self.add_main_section()

    def add_main_section(self):
        links = glob.glob("{}".format(self.pattern))
        names = [filename.rsplit("/", 1)[1].split(".html")[0] for filename in links]

        df = pd.DataFrame(
            {"names": names, "links": [link.split(os.sep, 1)[1] for link in links]}
        )
        df.sort_values(by="names")

        datatable = DataTable(df, "fastqc", index=False)
        datatable.datatable.set_links_to_column("links", "names")

        datatable.datatable.datatable_options = {
            "scrollX": "300px",
            "pageLength": 30,
            "scrollCollapse": "true",
            "dom": "rtpB",
            "paging": "false",
            "buttons": ["copy", "csv"],
        }
        js = datatable.create_javascript_function()
        html_tab = datatable.create_datatable()

        html = "{} {}".format(html_tab, js)

        self.sections.append(
            {
                "name": "FastQC report(s)",
                "anchor": "fastqc",
                "content": "<p> Here below are link(s) to original FastQC report. "
                "Please click on one of the links to jump to the main "
                "report.  {} </p>".format(html),
            }
        )
