# coding: utf-8
#
#  This file is part of Sequana software
#
#  Copyright (c) 2016 - Sequana Development Team
#
#  File author(s):
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

import pandas as pd

from sequana import bedtools
from sequana.modules_report.base_module import SequanaBaseModule
from sequana.utils import config
from sequana.utils.datatables_js import DataTable, DataTableFunction
from sequana.plots.canvasjs_linegraph import CanvasJSLineGraph
from sequana.modules_report.summary import SummaryModule

import colorlog

logger = colorlog.getLogger(__name__)


__all__ = ["CoverageModule", "ChromosomeCoverageModule"]


class CoverageModule(SequanaBaseModule):
    """Write HTML report of coverage analysis. This class takes either a
    :class:`GenomeCov` instances or a csv file where analysis are stored.
    """

    def __init__(self, data, region_window=200000):
        """.. rubric:: constructor

        :param data: it can be a csv filename created by sequana_coverage or a
            :class:`bedtools.GenomeCov` object.
        :param region_window:
        """
        super().__init__()
        self.region_window = region_window

        if isinstance(data, bedtools.GenomeCov):
            self.bed = data
        else:
            raise TypeError

        try:
            html_list = self.create_chromosome_reports()
        except TypeError:
            msg = (
                "Data must be either a csv file or a :class:`GenomeCov` "
                "instance where zscore is computed."
            )
            raise TypeError(msg)

        self.title = "Main coverage report ({0})".format(config.sample_name)
        self.intro = (
            "<p>Report the coverage of your sample ({0}) to check the "
            "quality of your mapping and to highlight regions of "
            "interest (under and over covered).</p>".format(config.sample_name)
        )
        self.create_report_content(html_list)
        self.create_html("sequana_coverage.html")

    def create_report_content(self, html_list):
        self.sections = list()
        self.create_chromosome_table(html_list)

    def create_chromosome_table(self, html_list):
        """Create table with links to chromosome reports"""
        df = pd.DataFrame(
            [
                [chrom.chrom_name, chrom.get_size(), chrom.DOC, chrom.CV, page]
                for chrom, page in zip(self.bed.chr_list, html_list)
            ],
            columns=["chromosome", "size", "mean_coverage", "coef_variation", "link"],
        )
        datatable = DataTable(df, "chrom")
        datatable.datatable.datatable_options = {
            "pageLength": 15,
            "dom": "Bfrtip",
            "buttons": ["copy", "csv"],
        }
        datatable.datatable.set_links_to_column("link", "chromosome")
        js = datatable.create_javascript_function()
        html_table = datatable.create_datatable(float_format="%.3g")
        self.sections.append(
            {
                "name": "Chromosomes",
                "anchor": "chromosomes",
                "content": "<p>Link to coverage analysis report for each chromosome. "
                "Size, mean coverage and coefficient of variation are reported"
                " in the table below.</p>\n{0}\n{1}".format(js, html_table),
            }
        )

    def create_chromosome_reports(self):
        """Create HTML report for each chromosome present in data."""
        # FIXME: why bed[0] (i.e. first chromosome)
        datatable_js = CoverageModule.init_roi_datatable(self.bed[0])
        chrom_output_dir = config.output_dir
        if not os.path.exists(chrom_output_dir):
            os.makedirs(chrom_output_dir)

        page_list = []
        for chrom in self.bed:
            logger.info("Creating coverage report {}".format(chrom.chrom_name))
            chrom_report = ChromosomeCoverageModule(
                chrom, datatable_js, region_window=self.region_window
            )
            page_list.append(chrom_report.html_page)
        return page_list

    # a static method because we need it in the coverage standalone
    # to initiate the datatables
    def init_roi_datatable(rois):
        """Initiate :class:`DataTableFunction` to create table to link each
        row with sub HTML report. All table will have the same appearance.
        We can therefore initialise the roi once for all.

        :param rois: can be a ROIs from ChromosomeCov instance or a simple
            dataframe
        """
        # computed
        try:
            df = rois.df.copy()
        except:
            df = rois.copy()
        df["link"] = ""
        # set datatable options
        datatable_js = DataTableFunction(df, "roi")
        if "start" in df.columns:
            datatable_js.set_links_to_column("link", "start")
        if "end" in df.columns:
            datatable_js.set_links_to_column("link", "end")
        datatable_js.datatable_options = {
            "scrollX": "true",
            "pageLength": 15,
            "scrollCollapse": "true",
            "dom": "Bfrtip",
            "buttons": ["copy", "csv"],
        }
        return datatable_js


class ChromosomeCoverageModule(SequanaBaseModule):
    """Write HTML report of coverage analysis for each chromosome. It is
    created by CoverageModule.
    """

    def __init__(
        self, chromosome, datatable, region_window=200000, options=None, command=""
    ):
        """

        :param chromosome:
        :param datatable:
        :param directory:
        :param int region_window: length of the sub coverage plot
        :param options: should contain "W", "k", "circular"

        """
        super().__init__()
        self.region_window = region_window

        directory = chromosome.chrom_name
        # to define where are css and js
        if directory in {None, "."}:
            self.path = ""
            directory = "."
        else:
            self.path = "../"

        self.chromosome = chromosome
        self.datatable = datatable
        self.command = command
        self.command += "\nSequana version: {}".format(config.version)
        self.title = "Coverage analysis of chromosome {0}".format(
            self.chromosome.chrom_name
        )

        self.intro = (
            "<p>The genome coverage analysis of the chromosome "
            "<b>{0}</b>.</p>".format(self.chromosome.chrom_name)
        )

        self.create_report_content(directory, options=options)
        self.html_page = "{0}{1}{2}.cov.html".format(
            directory, os.sep, self.chromosome.chrom_name
        )

        self.create_html(self.html_page)

    def create_report_content(self, directory, options=None):
        """Generate the sections list to fill the HTML report."""
        self.sections = list()

        if self.chromosome._mode == "memory":
            # nothing to do, all computations should be already available
            # and in memory
            rois = self.chromosome.get_rois()
        elif self.chromosome._mode == "chunks":
            # we need to reset the data to the first chunk
            # and compute the median and zscore. So, first, we save the entire
            # data set
            # self.chromosome.reset()
            # self.chromosome.running_median(options['W'], circular=options['circular'])
            # self.chromosome.compute_zscore(options['k'])
            # We mus set the ROI manually
            rois = options["ROIs"]

        if self.chromosome.DOC > 0:
            self.coverage_plot()
            if self.chromosome._mode == "memory":
                links = self.subcoverage(rois, directory)
            else:
                links = None
        else:
            links = None
        self.basic_stats()

        if self.chromosome.DOC:
            self.regions_of_interest(rois, links)
            self.coverage_barplot()
            if "gc" in self.chromosome.df.columns:
                self.gc_vs_coverage()
            self.normalized_coverage()
            self.zscore_distribution()
        self.add_command()

    def add_command(self):
        self.sections.append(
            {
                "name": "Command",
                "anchor": "command",
                "content": ("<p>Command used: <pre>{}</pre>.</p>".format(self.command)),
            }
        )

    def _add_large_data_section(self, name, anchor):
        self.sections.append(
            {
                "name": name,
                "anchor": anchor,
                "content": (
                    "<p>Large data sets (--chunk-size argument "
                    "used), skipped plotting.</p>"
                ),
            }
        )

    def coverage_plot(self):
        """Coverage section."""
        if self.chromosome._mode == "chunks":
            self._add_large_data_section("Coverage", "coverage")
            return

        image = self.create_embedded_png(
            self.chromosome.plot_coverage, input_arg="filename"
        )

        self.sections.append(
            {
                "name": "Coverage",
                "anchor": "coverage",
                "content": "<p>The following figures shows the per-base coverage along the"
                " reference genome (black line). The blue line indicates the "
                "running median. From the normalised coverage, we estimate "
                "z-scores on a per-base level. The red lines indicates the "
                "z-scores at plus or minus N standard deviations, where N is "
                "chosen by the user. (default:4). Only a million point are "
                "shown. This may explain some visual discrepancies with. </p>\n{0}".format(
                    image
                ),
            }
        )

    def coverage_barplot(self):
        """Coverage barplots section."""
        if self.chromosome._mode == "chunks":
            self._add_large_data_section("Coverage histogram", "cov_barplot")
            return

        image1 = self.create_embedded_png(
            self.chromosome.plot_hist_coverage, input_arg="filename", style="width:45%"
        )
        image2 = self.create_embedded_png(
            self.chromosome.plot_hist_coverage,
            input_arg="filename",
            style="width:45%",
            logx=False,
        )
        self.sections.append(
            {
                "name": "Coverage histogram",
                "anchor": "cov_barplot",
                "content": "<p>The following figures contain the histogram of the genome "
                "coverage. The X and Y axis being in log scale in the left panel"
                "while only the Y axis is in log scale in the right panel.</p>\n"
                "{0}\n{1}".format(image1, image2),
            }
        )

    def subcoverage(self, rois, directory):
        """Create subcoverage reports to have access to a zoomable line plot.

        :params rois:
        :param directory: directory name for the chromsome

        This method create sub reports for each region of 200,000 bases (can be
        changed). Usually, it starts at position 0 so reports will be stored
        in e.g. for a genome of 2,300,000 bases::

            chromosome_name/chromosome_name_0_200000.html
            chromosome_name/chromosome_name_200000_400000.html
            ...
            ...
            chromosome_name/chromosome_name_2000000_2200000.html
            chromosome_name/chromosome_name_2200000_2300000.html

        Note that if the BED file positions does not start at zero, then
        names will take care of that.

        """
        # an aliases
        W = self.region_window
        name = self.chromosome.chrom_name
        chrom = self.chromosome
        N = len(self.chromosome)

        # position does not always start at position zero
        shift = self.chromosome.df.pos.iloc[0]
        maxpos = self.chromosome.df.pos.iloc[-1]

        # create directory

        chrom_output_dir = os.sep.join([config.output_dir, str(directory), str(name)])
        if not os.path.exists(chrom_output_dir):
            os.makedirs(chrom_output_dir)

        # create the combobox to link toward different sub coverage
        # Here, we should (1) get the length of the data and (2)
        # figure out the boundary. Indeed, you can imagine a BED file
        # that does not start at position zero, but from POS1>0 to POS2

        links = [
            "{0}/{0}_{1}_{2}.html".format(name, i, min(i + W, maxpos))
            for i in range(shift, shift + N, W)
        ]
        intra_links = (
            "{0}_{1}_{2}.html".format(name, i, min(i + W, maxpos))
            for i in range(shift, shift + N, W)
        )

        combobox = self.create_combobox(links, "sub", True)
        combobox_intra = self.create_combobox(intra_links, "sub", False)
        datatable = self._init_datatable_function(rois)

        # break the chromosome as pieces of 200,000 bp
        for i in range(shift, shift + N, W):
            SubCoverageModule(
                chrom, rois, combobox_intra, datatable, i, min(i + W, maxpos), directory
            )

        self.sections.append(
            {"name": "Subcoverage", "anchor": "subcoverage", "content": combobox}
        )
        return links

    def _init_datatable_function(self, rois):
        """Initiate :class:`DataTableFunction` to create table to link each
        row with sub HTML report. All table will have the same appearance. So,
        let's initiate its only once.
        """
        datatable_js = DataTableFunction(rois.df, "roi")
        datatable_js.datatable_options = {
            "scrollX": "true",
            "pageLength": 15,
            "scrollCollapse": "true",
            "dom": "Bfrtip",
            "buttons": ["copy", "csv"],
        }
        return datatable_js

    def basic_stats(self):
        """Basics statistics section."""
        li = "<li><b>{0}</b> ({1}): {2:.2f}</li>"
        stats = self.chromosome.get_stats()

        description = {
            "BOC": "breadth of coverage: the proportion (in %s) "
            + " of the genome covered by at least one read.",
            "CV": "the coefficient of variation.",
            "DOC": "the sequencing depth (Depth of Coverage), that is the "
            + "average the genome coverage.",
            "MAD": "median of the absolute median deviation defined as median(|X-median(X)|).",
            "Median": "Median of the coverage.",
            "STD": "standard deviation.",
            "GC": "GC content (%)",
        }

        data = [
            li.format(key, description[key], stats[key])
            for key in ["BOC", "CV", "DOC", "MAD", "Median", "STD", "GC"]
            if key in stats.keys()
        ]

        stats = "<ul>{0}</ul>".format("\n".join(data))
        self.sections.append(
            {
                "name": "Basic stats",
                "anchor": "basic_stats",
                "content": "<p>Here are some basic statistics about the "
                "genome coverage.</p>\n{0}".format(stats),
            }
        )

    def regions_of_interest(self, rois, links):
        """Region of interest section."""

        def connect_link(x):
            for link in links:
                _, x1, x2 = link.rsplit(os.sep)[1].rstrip(".html").rsplit("_", 2)
                x1 = int(x1)
                x2 = int(x2)
                if x >= x1 and x <= x2:
                    return link
            # for the case where the data is fully stored in memory, we must
            # find all events !
            if self.chromosome._mode == "memory" and self.chromosome.binning == 1:
                raise Exception("{} position not in the range of reports".format(x))

        if links:
            links_list = [connect_link(n) for n in rois.df["start"]]
        else:
            links_list = [None for n in rois.df["start"]]

        rois.df["link"] = links_list
        # create datatable
        low_roi = rois.get_low_rois()
        high_roi = rois.get_high_rois()

        datatable = CoverageModule.init_roi_datatable(low_roi)

        datatable.set_links_to_column("link", "chr")

        js = datatable.create_javascript_function()
        lroi = DataTable(low_roi, "lroi", datatable)
        hroi = DataTable(high_roi, "hroi", datatable)
        html_low_roi = lroi.create_datatable(float_format="%.3g")
        html_high_roi = hroi.create_datatable(float_format="%.3g")
        rois.df.drop("link", 1, inplace=True)
        roi_paragraph = (
            "<p>Regions with a z-score {0}er than {1:.2f} and at "
            "least one base with a z-score {0}er than {2:.2f} are detected."
            "There are {3} {0} regions of interest."
            "</p>"
        )
        low_paragraph = roi_paragraph.format(
            "low",
            self.chromosome.thresholds.low2,
            self.chromosome.thresholds.low,
            len(low_roi),
        )
        high_paragraph = roi_paragraph.format(
            "high",
            self.chromosome.thresholds.high2,
            self.chromosome.thresholds.high,
            len(high_roi),
        )

        self.sections.append(
            {
                "name": "Regions Of Interest (ROI)",
                "anchor": "roi",
                "content": "{4}\n"
                "<p>The following tables give regions of "
                "interest detected by sequana. Here are the definitions of the "
                "columns:</p>\n"
                "<ul><li>mean_cov: the average of coverage</li>\n"
                "<li>mean_rm: the average of running median</li>\n"
                "<li>mean_zscore: the average of zscore</li>\n"
                "<li>max_zscore: the higher zscore contains in the region</li>"
                "</ul>\n"
                "<h3>Low coverage region</h3>\n{0}\n{1}\n"
                "<h3>High coverage region</h3>\n{2}\n{3}\n".format(
                    low_paragraph, html_low_roi, high_paragraph, html_high_roi, js
                ),
            }
        )

    def gc_vs_coverage(self):
        """3 dimensional plot of GC content versus coverage."""
        image = self.create_embedded_png(
            self.chromosome.plot_gc_vs_coverage, input_arg="filename"
        )
        corr = self.chromosome.get_gc_correlation()
        self.sections.append(
            {
                "name": "Coverage vs GC content",
                "anchor": "cov_vs_gc",
                "content": "<p>The correlation coefficient between the coverage and GC "
                "content is <b>{0:.3g}</b> with a window size of {1}bp.</p>\n"
                "{2}\n"
                "<p>Note: the correlation coefficient has to be between -1.0 "
                "and 1.0. A coefficient of 0 means no correlation, while a "
                "coefficient of -1 or 1 means an existing "
                "correlation between GC and Coverage</p>".format(
                    corr, self.chromosome.bed.gc_window_size, image
                ),
            }
        )

    def normalized_coverage(self):
        """Barplot of normalized coverage section."""
        if self.chromosome._mode == "chunks":
            self._add_large_data_section("Normalised coverage", "normalised_coverage")
            return

        image = self.create_embedded_png(
            self.chromosome.plot_hist_normalized_coverage, input_arg="filename"
        )
        self.sections.append(
            {
                "name": "Normalised coverage",
                "anchor": "normalised_coverage",
                "content": "<p>Distribution of the normalised coverage with predicted "
                "Gaussian. The red line should be followed the trend of the "
                "barplot.</p>\n{0}".format(image),
            }
        )

    def zscore_distribution(self):
        """Barplot of zscore distribution section."""
        if self.chromosome._mode == "chunks":
            self._add_large_data_section("Z-Score distribution", "zscore_distribution")
            return

        image = self.create_embedded_png(
            self.chromosome.plot_hist_zscore, input_arg="filename"
        )
        self.sections.append(
            {
                "name": "Z-Score distribution",
                "anchor": "zscore_distribution",
                "content": "<p>Distribution of the z-score (normalised coverage); You "
                "should see a Gaussian distribution centered around 0. The "
                "estimated parameters are mu={0:.2f} and sigma={1:.2f}.</p>\n"
                "{2}".format(
                    self.chromosome.best_gaussian["mu"],
                    self.chromosome.best_gaussian["sigma"],
                    image,
                ),
            }
        )


class SubCoverageModule(SequanaBaseModule):
    """Write HTML report of subsection of chromosome with a javascript
    coverage plot.
    """

    def __init__(self, chromosome, rois, combobox, datatable, start, stop, directory):
        super().__init__()

        if directory == ".":
            self.path = "../"
        else:
            self.path = "../../"

        self.chromosome = chromosome
        self.rois = rois
        self.combobox = combobox
        self.datatable = datatable

        self.start = start
        self.stop = stop
        self.title = (
            "Coverage analysis of chromosome {0}<br>"
            "positions {1} and {2}".format(self.chromosome.chrom_name, start, stop)
        )
        self.create_report_content()
        self.create_html(
            "{0}{4}{1}{4}{1}_{2}_{3}.html".format(
                directory, self.chromosome.chrom_name, start, stop, os.sep
            )
        )

    def create_report_content(self):
        """Generate the sections list to fill the HTML report."""
        self.sections = list()

        self.canvasjs_line_plot()
        self.regions_of_interest()

    def canvasjs_line_plot(self):
        """Create the CanvasJS line plot section."""
        # set column of interest and create the csv
        x_col = "pos"
        y_col = ("cov", "mapq0", "gc")
        columns = self.chromosome.df.columns
        y_col = [n for n in y_col if n in columns]
        csv = self.chromosome.to_csv(
            start=self.start,
            stop=self.stop,
            columns=[x_col] + y_col,
            index=False,
            float_format="%.3g",
        )

        # create CanvasJS stuff
        cjs = CanvasJSLineGraph(csv, "cov", x_col, y_col)
        # set options
        cjs.set_options(
            {"zoomEnabled": "true", "zoomType": "x", "exportEnabled": "true"}
        )
        # set title
        cjs.set_title("Genome Coverage")
        # set legend
        cjs.set_legend(
            {
                "verticalAlign": "bottom",
                "horizontalAlign": "center",
                "cursor": "pointer",
            },
            hide_on_click=True,
        )
        # set axis
        cjs.set_axis_x(
            {
                "title": "Position (bp)",
                "labelAngle": 30,
                "minimum": self.start,
                "maximum": self.stop,
            }
        )
        cjs.set_axis_y({"title": "Coverage (Count)"})
        cjs.set_axis_y2(
            {
                "title": "GC content (ratio)",
                "minimum": 0,
                "maximum": 1,
                "lineColor": "#FFC425",
                "titleFontColor": "#FFC425",
                "labelFontColor": "#FFC425",
            }
        )
        # set datas
        cjs.set_data(
            index=0,
            data_dict={
                "type": "line",
                "name": "Coverage",
                "showInLegend": "true",
                "color": "#5BC0DE",
                "lineColor": "#5BC0DE",
            },
        )
        try:
            i = y_col.index("mapq0")
            cjs.set_data(
                index=i,
                data_dict={
                    "type": "line",
                    "name": "Filtered coverage",
                    "showInLegend": "true",
                    "color": "#D9534F",
                    "lineColor": "#D9534F",
                },
            )
        except ValueError:
            pass
        try:
            i = y_col.index("gc")
            cjs.set_data(
                index=i,
                data_dict={
                    "type": "line",
                    "axisYType": "secondary",
                    "name": "GC content",
                    "showInLegend": "true",
                    "color": "#FFC425",
                    "lineColor": "#FFC425",
                },
            )
        except ValueError:
            pass
        # create canvasJS
        html_cjs = cjs.create_canvasjs()
        self.sections.append(
            {
                "name": "Interactive coverage plot",
                "anchor": "iplot",
                "content": ("{0}{1}\n".format(self.combobox, html_cjs)),
            }
        )

    def regions_of_interest(self):
        """Region of interest section."""
        subseq = [self.start, self.stop]

        low_roi = self.rois.get_low_rois(subseq)
        high_roi = self.rois.get_high_rois(subseq)
        js = self.datatable.create_javascript_function()
        lroi = DataTable(low_roi, "lroi", self.datatable)
        hroi = DataTable(high_roi, "hroi", self.datatable)
        html_low_roi = lroi.create_datatable(float_format="%.3g")
        html_high_roi = hroi.create_datatable(float_format="%.3g")
        roi_paragraph = (
            "<p>Regions with a z-score {0}er than {1:.2f} and at "
            "least one base with a z-score {0}er than {2:.2f} are detected as "
            "{0} coverage region. Thus, there are {3} {0} coverage regions "
            "between the position {4} and the position {5}</p>"
        )
        low_paragraph = roi_paragraph.format(
            "low",
            self.chromosome.thresholds.low2,
            self.chromosome.thresholds.low,
            len(low_roi),
            self.start,
            self.stop,
        )
        high_paragraph = roi_paragraph.format(
            "high",
            self.chromosome.thresholds.high2,
            self.chromosome.thresholds.high,
            len(high_roi),
            self.start,
            self.stop,
        )

        self.sections.append(
            {
                "name": "Regions Of Interest (ROI)",
                "anchor": "roi",
                "content": "{4}\n"
                "<p>Running median is the median computed along the genome "
                "using a sliding window. The following tables give regions of "
                "interest detected by sequana. Here are some definition of the "
                "table's columns:</p>\n"
                "<ul><li><b>mean_cov</b>: the average of coverage</li>\n"
                "<li><b>mean_rm</b>: the average of running median</li>\n"
                "<li><b>mean_zscore</b>: the average of zscore</li>\n"
                "<li><b>max_zscore</b>: the higher zscore contains in the "
                "region</li>\n"
                "<li><b>log2_ratio</b>:log2(mean_cov/mean_rm)</li></ul>\n"
                "<h3>Low coverage region</h3>\n{0}\n{1}\n"
                "<h3>High coverage region</h3>\n{2}\n{3}\n".format(
                    low_paragraph, html_low_roi, high_paragraph, html_high_roi, js
                ),
            }
        )
