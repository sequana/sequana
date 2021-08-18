# -*- coding: utf-8 -*-
#
#  This file is part of Sequana software
#
#  Copyright (c) 2016 - Sequana Development Team
#
#  File author(s):
#      Thomas Cokelaer <thomas.cokelaer@pasteur.fr>
#      Dimitri Desvillechabrol <dimitri.desvillechabrol@pasteur.fr>,
#          <d.desvillechabrol@gmail.com>
#      Rachel Legendre <rachel.legendre@pasteur.fr>
#
#  Distributed under the terms of the 3-clause BSD license.
#  The full license is in the LICENSE file, distributed with this software.
#
#  website: https://github.com/sequana/sequana
#  documentation: http://sequana.readthedocs.io
#
##############################################################################
"""Report dedicated to BAM file

.. autosummary::

    TRFModule

"""
import os

from sequana.lazy import pandas as pd
from sequana.modules_report.base_module import SequanaBaseModule

from sequana import TRF

from sequana.lazy import pylab

from sequana.utils.datatables_js import DataTable
from sequana.utils.df2html import df2html

__all__ = ["TRFModule"]


class TRFModule(SequanaBaseModule):
    """Report dedicated to TRF file"""

    def __init__(self, trf_input, output_filename=None):
        super().__init__()

        self.trf = TRF(trf_input)
        self.title = "TRF Report"
        self.create_report_content()
        self.create_html(output_filename)

    def create_report_content(self):
        self.sections = list()
        self.add_info()
        self.add_table()
        self.add_images_section()

    def add_info(self):

        html = "<br>".join(self.trf.__repr__().split("\n"))

        self.sections.append({"name": "Information", "anchor": "info", "content": html})

    def add_table(self):

        datatable = DataTable(self.trf.df, "result", index=True)
        datatable.datatable.datatable_options = {
            "scrollX": "300px",
            "pageLength": 15,
            "scrollCollapse": "true",
            "dom": "tBifp",
            "paging": "true",
            "buttons": ["copy", "csv"],
        }
        js = datatable.create_javascript_function()
        html_tab = datatable.create_datatable(float_format="%.3g")

        html = ""
        html += "{} {}".format(html_tab, js)

        self.sections.append(
            {"name": "TRF results", "anchor": "results", "content": html}
        )

    def add_images_section(self):
        style = "width:65%"
        import pylab

        pylab.ioff()

        def plotter1(filename):
            pylab.clf()
            self.trf.hist_entropy()
            pylab.savefig(filename)

        html1 = self.create_embedded_png(plotter1, "filename", style=style)

        def plotter2(filename):
            pylab.clf()
            self.trf.hist_period_size()
            pylab.savefig(filename)

        html2 = self.create_embedded_png(plotter2, "filename", style=style)

        self.sections.append(
            {"name": "Image", "anchor": "table", "content": html1 + html2}
        )
