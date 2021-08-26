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

    BAMQCModule

"""
import os

from sequana.lazy import pandas as pd
from sequana.modules_report.base_module import SequanaBaseModule

from sequana.bamtools import SAMFlags
from sequana import BAM

from sequana.lazy import pylab

from sequana.utils.datatables_js import DataTable

__all__ = ["BAMQCModule"]


class BAMQCModule(SequanaBaseModule):
    """Report dedicated to BAM file

    ::

        from sequana import sequana_data
        from sequana.modules_report.bamqc import BAMQCModule
        filename = sequana_data("test.bam")

        r = BAMQCModule(filename)
        r.create_html("test.html")

        # report/bam.html is now available

    .. todo:: right now, the computation is performed in the class. Ideally,
        we would like the computation to happen elsewhere, where a json is stored.
        The json would be the input to this class.
    """

    def __init__(self, bam_input, output_filename=None):
        super().__init__()

        self.bam_input = bam_input
        self.title = "Bam Report"
        self.create_report_content()
        self.create_html(output_filename)

    def create_report_content(self):
        self.sections = list()
        self.add_flag_section()
        self.add_images_section()

    def _computation(self):
        self.bam = BAM(self.bam_input)

        results = {}
        results["alignment_count"] = len(self.bam)

        # first, we store the flags
        df = self.bam.get_flags_as_df().sum()
        df = df.to_frame()
        df.columns = ["counter"]
        sf = SAMFlags()
        df["meaning"] = sf.get_meaning()
        df = df[["meaning", "counter"]]
        results["flags"] = df

        return results

        self.bam.plot_bar_flags(
            logy=False, filename=self.directory + os.sep + "bar_flags.png"
        )
        self.bam.plot_bar_mapq(filename=self.directory + os.sep + "bar_mapq.png")

    def add_flag_section(self):
        data = self._computation()
        df = data["flags"]

        datatable = DataTable(df, "flags", index=True)
        datatable.datatable.datatable_options = {
            "scrollX": "300px",
            "pageLength": 15,
            "scrollCollapse": "true",
            "dom": "tB",
            "paging": "false",
            "buttons": ["copy", "csv"],
        }
        js = datatable.create_javascript_function()
        html_tab = datatable.create_datatable(float_format="%.3g")

        html = ""
        html += "{} {}".format(html_tab, js)

        self.sections.append(
            {"name": "Flags information", "anchor": "flags", "content": html}
        )

    def add_images_section(self):
        style = "width:65%"
        import pylab

        pylab.ioff()

        def plotter1(filename):
            self.bam.plot_bar_flags(logy=True, filename=filename)

        html1 = self.create_embedded_png(plotter1, "filename", style=style)

        def plotter2(filename):
            self.bam.plot_bar_flags(logy=False, filename=filename)

        html2 = self.create_embedded_png(plotter2, "filename", style=style)

        def plotter3(filename):
            self.bam.plot_bar_mapq(filename=filename)

        html3 = self.create_embedded_png(plotter3, "filename", style=style)

        self.sections.append(
            {"name": "Image", "anchor": "table", "content": html1 + html2 + html3}
        )
