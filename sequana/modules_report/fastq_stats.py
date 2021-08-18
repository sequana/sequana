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


class FastQStatsModule(SequanaBaseModule):
    """Write HTML report of fastq stats analysis."""

    def __init__(
        self, input_directory, path_to_fastqc, output_filename=None, tag_R1="_R1_"
    ):
        """
        :param input_directory: where to find the json and boxplot image. The
            path where to find the data does not matter since the JSON and PNG
            will be embedded.
        :param path_to_fastqc: This must be provided by the user. This is the
            directory where will be found the original FastQC reports. This can
            be infered but is prone to error so for now, we must provide this
            argument.
        :param output_filename: if not provided, the HTML is not created.

        ::

            from sequana.modules_report.fastq_stats import FastQStatsModule
            ff = FastQStatsModule("./SAMPLE/fastq_stats_samples", "fastqc_samples",
                "test.html")

        """
        super().__init__()
        self.path_to_fastqc = path_to_fastqc
        self.directory = input_directory
        self.create_report_content()
        if output_filename:
            self.create_html(output_filename)

    def create_report_content(self):
        """Generate the sections list to fill the HTML report."""
        self.sections = list()
        self.add_stats()

    def _get_files(self, pattern):
        # !! need to sort the files so that R1 appears before R2
        filenames = sorted(glob.glob(self.directory + os.sep + pattern))
        if len(filenames) == 2:
            mode = "pe"
        elif len(filenames) == 1:
            mode = "se"
        elif len(filenames) == 0:
            return
        else:
            logger.warning(
                "FastQStatsModule: more than 2 files "
                "matched the pattern %s" % pattern
            )
            return
        return filenames, mode

    def get_stats(self):
        import pandas as pd

        filenames, mode = self._get_files("*.json")
        if mode == "pe":
            df1 = pd.read_json(filenames[0])
            df2 = pd.read_json(filenames[1])
            df = pd.concat([df1, df2])
            # Should have been sorted !
            df.index = ["R1", "R2"]
        else:
            df = pd.read_json(filenames[0])
            df.index = ["R1"]
        df = df[
            [
                "A",
                "C",
                "G",
                "T",
                "N",
                "n_reads",
                "mean quality",
                "GC content",
                "average read length",
                "total bases",
            ]
        ]
        for this in "ACGTN":
            df[this] /= df["total bases"]
            df[this] *= 100
        return df

    def _get_stats_section(self, tablename="stats"):
        self.df_stats = self.get_stats()
        filenames, mode = self._get_files("*boxplot.png")

        datatable = DataTable(self.df_stats, tablename, index=True)
        datatable.datatable.datatable_options = {
            "scrollX": "300px",
            "pageLength": 15,
            "scrollCollapse": "true",
            "dom": "rtpB",
            "paging": "false",
            "buttons": ["copy", "csv"],
        }
        js = datatable.create_javascript_function()
        html_tab = datatable.create_datatable(float_format="%.3g")

        html = """<p>The following table gives some basic statistics about the data before any filtering.
   The A, C, G, T, N columns report the percentage of each bases in the overall sequences.
   The GC content is provided in percentage as well. </p>
   <div>{} {}</div>
   <div>""".format(
            html_tab, js
        )

        html += """
   <p>The following figure(s) gives the average quality (red line) of raw reads
   (500,000 at max). The x-axis being the length of the reads. The yellow
   enveloppe gives the variation of the quality (1 standard deviation).</p>
   <p> Click on the image to jump to a full FastQC report.</p>"""

        if len(filenames) == 2:
            width = "49"
        else:
            width = "65"

        filename = os.path.split(filenames[0])[1].replace(
            "_boxplot.png", "_fastqc.html"
        )
        href = self.path_to_fastqc + os.sep + filename
        html += """
   <figure style="float:left; width:{}%; padding:0px; margin:0px;">
       <a href="{}">{}</a>
   <figcaption style="font-style:italic">Fig1: R1 reads</figcaption>
   </figure>""".format(
            width, href, self.png_to_embedded_png(filenames[0])
        )

        if len(filenames) == 2:
            filename = os.path.split(filenames[1])[1].replace(
                "_boxplot.png", "_fastqc.html"
            )
            href = self.path_to_fastqc + os.sep + filename
            html += """
   <figure style="float:right; width:{}%; padding:0px; margin:0px;">
       <a href="{}">{}</a>
   <figcaption style="font-style:italic">Fig2: R2 reads</figcaption>
   </figure>""".format(
                width, href, self.png_to_embedded_png(filenames[1])
            )

        return html

    def add_stats(self):
        html = self._get_stats_section()
        self.sections.append(
            {"name": "Stats inputs", "anchor": "stats", "content": html}
        )
