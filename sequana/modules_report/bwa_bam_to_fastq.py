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

import colorlog

logger = colorlog.getLogger(__name__)


from sequana.utils.datatables_js import DataTable


class BWABAMtoFastQModule(SequanaBaseModule):
    """Write HTML report of BWA mapping (phix)"""

    def __init__(self, input_directory, output_filename=None):
        """
        :param input_directory: the directory of the bwa_bam_to_fastq output
        :param output_filename: if not provided, the HTML is not created.

        """
        super().__init__()
        self.directory = input_directory + os.sep
        self.create_report_content()
        if output_filename:
            self.create_html(output_filename)

    def create_report_content(self):
        """Generate the sections list to fill the HTML report."""
        self.sections = list()
        self.add_stats()

    def _get_html_stats(self):
        from sequana.tools import StatsBAM2Mapped
        from easydev import precision

        data = StatsBAM2Mapped(self.directory + "bwa_mem_stats.json").data
        html = "Reads with Phix: %s %%<br>" % precision(data["contamination"], 3)

        # add HTML table
        if "R2_mapped" in data.keys():
            df = pd.DataFrame(
                {
                    "R1": [data["R1_mapped"], data["R1_unmapped"]],
                    "R2": [data["R2_mapped"], data["R2_unmapped"]],
                }
            )
        else:
            df = pd.DataFrame({"R1": [data["R1_mapped"], data["R1_unmapped"]]})
        df.index = ["mapped", "unmapped"]

        datatable = DataTable(df, "bwa_bam")
        datatable.datatable.datatable_options = {
            "scrollX": "300px",
            "pageLength": 30,
            "scrollCollapse": "true",
            "dom": "irtpB",
            "paging": "false",
            "buttons": ["copy", "csv"],
        }
        js = datatable.create_javascript_function()
        html_tab = datatable.create_datatable(float_format="%.3g")
        # html += "{} {}".format(html_tab, js)

        html += "Unpaired: %s <br>" % data["unpaired"]
        html += "duplicated: %s <br>" % data["duplicated"]
        return html

    def _get_html_mapped_stats(self):
        html = ""
        return html

    def add_stats(self):
        html1 = self._get_html_stats()
        html2 = self._get_html_mapped_stats()
        self.sections.append(
            {"name": "Stats inputs", "anchor": "stats", "content": html1 + html2}
        )
