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
import json

from sequana.modules_report.base_module import SequanaBaseModule
from sequana.utils import config

from sequana.lazy import pandas as pd
from sequana.lazy import pylab
from sequana.lazy import pandas

from sequana.utils.datatables_js import DataTable


class PacbioInputBAMModule(SequanaBaseModule):
    """Write HTML report of Pacbio input bam.

    Input summary JSON file must contains these links:

    - images/hist_read_length
    - images/hist_gc_content
    - images/hist_zmw

    to PNG files and the **stats** dictionary created with
    :meth:`sequana.pacbio.BAMPacbio.stats

    """

    def __init__(self, summary, output_filename=None):
        """
        :param input:
        """
        super().__init__()
        # Read the data
        self.title = "Pacbio QC (input BAM)"
        self.name = "sequana_summary_pacbio_qc"
        with open(summary, "r") as fh:
            self.summary = json.load(fh)

        assert "name" in self.summary.keys() and self.summary["name"] == self.name

        self.create_report_content()
        if output_filename:
            self.create_html(output_filename)

    def create_report_content(self):
        """Generate the sections list to fill the HTML report."""
        self.sections = list()

        self.add_stats()

        # we use a list to keep this order
        self.key_images = [
            "hist_read_length",
            "hist_gc_content",
            "gc_vs_length",
            "hist_snr",
            "hist_zmw",
        ]
        for this in self.key_images:
            self.add_png(this)

    def add_stats(self):
        df = pd.Series(self.summary["read_stats"]).to_frame().T
        df.index = ["read length stats"]
        table = DataTable(df, "table", index=True)
        table.datatable.datatable_options = {
            "scrollX": "300px",
            "pageLength": 30,
            "scrollCollapse": "true",
            "dom": "t",
            "paging": "false",
            "buttons": ["copy", "csv"],
        }
        js = table.create_javascript_function()
        # IMPORTANT: here conversion to integer with %d
        # to round and make integer. !! The GC is therefore
        # converted to integer as well.
        html_tab = table.create_datatable(float_format="%d")
        html = "{} {}".format(html_tab, js)

        self.sections.append(
            {"name": "Basic stats on read length", "anchor": "table", "content": html}
        )

    def add_png(self, key):
        text = ""
        if key == "hist_read_length":
            title = "Histogram read length"
        elif key == "hist_gc_content":
            title = "Histogram GC content"
        elif key == "hist_zmw":
            title = "Histogram ZMW"
        elif key == "gc_vs_length":
            title = "GC vs length"
        elif key == "hist_snr":
            title = "Histogram SNR"
            text = (
                "SNR histogram of A, C, G, T bases from the raw data if "
                "present in the BAM file, otherwise you will see a 'no "
                "data' icon"
            )
        else:
            return

        html = self.png_to_embedded_png(self.summary["images"][key])
        html = """<figure style="float:center; width:89%; padding:0px; margin:0px;">
        {}
    <figcaption style="font-style:italic"></figcaption>
    </figure>""".format(
            html, html
        )

        self.sections.append(
            {
                "name": title,
                "anchor": key,
                "content": "<p>{}</p><br>".format(text) + html,
            }
        )
