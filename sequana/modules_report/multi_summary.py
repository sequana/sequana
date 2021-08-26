#
#  Copyright (c) 2016-2021 - Sequana Development Team
#
#  The full license is in the LICENSE file, distributed with this software.
#
#  website: https://github.com/sequana/sequana
#  documentation: http://sequana.readthedocs.io
#
##############################################################################
import glob
import json

from sequana.modules_report.base_module import SequanaBaseModule

from easydev import DevTools
from sequana.resources.canvas.bar import CanvasBar

from sequana.lazy import pandas as pd

import colorlog

logger = colorlog.getLogger(__name__)


__all__ = ["ReadSummary", "MultiSummary"]


class ReadSummary(object):
    def __init__(self, filename):
        self.filename = filename
        self.data = json.load(open(self.filename, "r"))

    def get_phix_percent(self):
        return self.data["phix_section"]["contamination"]

    def get_cutadapt_stats(self):
        return self.data["cutadapt_json"]

    def get_fastq_stats_samples(self):
        return self.data["fastq_stats_samples_json"]

    def get_mean_quality_samples(self):
        this = self.data["fastq_stats_samples_json"]
        return this["mean quality"]["R1"]

    def get_nreads_raw(self):
        this = self.data["fastq_stats_samples_json"]
        return this["n_reads"]["R1"]

    def get_average_read_length(self):
        this = self.data["fastq_stats_samples_json"]
        return this["average read length"]["R1"]

    def get_trimming_percent(self):
        d = self.get_cutadapt_stats()
        try:
            trimming = d["percent"]["Pairs too short"]
        except KeyError:
            trimming = d["percent"]["Reads too short"]
        return float(
            trimming.replace("%", "").replace("(", "").replace(")", "").strip()
        )

    def get_adapters_percent(self):
        d = self.get_cutadapt_stats()
        try:
            read1 = d["percent"]["Read1 with adapters"]
            # e.g. (8.3%)
            read1 = float(read1.strip("%)").strip("("))
            try:
                read2 = d["percent"]["Read2 with adapters"]
                read2 = float(read2.strip("%)").strip("("))
                read1 = (read1 + read2) / 2.0  # FIXME crude approximation
            except:
                pass  # single-end
        except:
            read1 = d["percent"]["Reads with adapters"]
            read1 = float(read1.strip("%)").strip("("))
        return read1

    def get_read1_with_adapters_percent(self):
        return self._get_read_with_adapters_percent("1")

    def get_read2_with_adapters_percent(self):
        return self._get_read_with_adapters_percent("2")

    def _get_read_with_adapters_percent(self, tag):
        d = self.get_cutadapt_stats()
        trimming = d["percent"]["Read%s with adapters" % tag]
        trimming = trimming.strip()
        for this in [",", "(", ")", "%"]:
            trimming = trimming.replace(this, "")
        trimming = float(trimming)
        return trimming

    def get_output_total_reads(self):
        d = self.get_cutadapt_stats()
        try:
            trimming = d["Number of reads"]["Pairs kept"]
        except:
            trimming = d["Number of reads"]["Reads kept"]
        try:  # previous version stored strings in the json; TODO add test
            trimming = trimming.strip()
            for this in [",", "(", ")", "%"]:
                trimming = trimming.replace(this, "")
            trimming = int(trimming)
        except:
            pass
        return trimming


class MultiSummary(SequanaBaseModule):
    """
    Used by the pipelines to create a summary based on the content of the
    directory. Also used by the standalone application, in which case
    config and pipeline files are not required.

    For developers:

    1. In class Summary:

    2 In class SequanaMultipleSummary:

        try:self.populate_gc()
        except:pass
        def populate_gc():
            do something
        def get_gc(self):
            return [x.get_gc() for x in self.summaries]

    3: update the jinja file report_multiple_summay

    """

    def __init__(
        self, pattern="**/summary.json", output_filename=None, verbose=True, **kargs
    ):
        super().__init__()

        from sequana import logger

        logger.setLevel("INFO")
        if verbose is False:
            logger.setLevel("WARNING")

        logger.info(
            "Sequana Quality control Summary will not be maintained. Please use sequana_fastqc for QCs, sequana_multitax for Taxonomy."
        )
        self.title = "Sequana multiple summary"
        self.devtools = DevTools()

        self.filenames = list(glob.iglob(pattern, recursive=True))

        self.summaries = [ReadSummary(filename) for filename in self.filenames]
        self.projects = [s.data["project"] for s in self.summaries]

        self.create_report_content()
        self.create_html(output_filename)

    def create_report_content(self):
        self.sections = list()
        self.add_section()

    def add_section(self):
        logger.info("Found %s projects/samples/ directories" % len(self.summaries))
        for filename in self.filenames:
            logger.info(filename)

        self.jinja = {}

        self.jinja[
            "canvas"
        ] = '<script type="text/javascript" src="js/canvasjs.min.js"></script>'
        self.jinja[
            "canvas"
        ] += """<script type="text/javascript">
            window.onload = function () {"""

        # Information to put on top of the page (added later in a module.intro)
        # We should get the link name from the project name contained in the json
        links = [
            {"href": filename.replace(".json", ".html"), "caption": project}
            for filename, project in zip(self.filenames, self.projects)
        ]
        introhtml = "<div><b>Number of samples:</b>{}</div>".format(len(self.summaries))

        self.jinja["sections"] = []

        # This will used to stored all information
        self.df = {}

        # The order does not matter here, everything is done in JINJA
        try:
            self.populate_nreads_raw()
        except Exception as err:
            print(err)

        try:
            self.populate_phix()
        except Exception as err:
            logger.warning("multi_summary: skip phix")

        try:
            self.populate_gc_samples()
        except Exception as err:
            logger.debug("multi_summary: skip gc samples")

        try:
            self.populate_trimming()
        except Exception as err:
            logger.debug("multi_summary: skip trimming")

        try:
            self.populate_mean_quality()
        except Exception as err:
            logger.debug("multi_summary: skip mean quality")

        try:
            self.populate_adapters()
        except Exception as err:
            logger.debug("multi_summary: skip adapters")

        try:
            self.populate_output_total_reads()
        except Exception as err:
            logger.debug("multi_summary: skip total reads")

        # Now we have all data in df as dictionaries. Let us merge them together

        keys = list(self.df.keys())
        if len(keys) >= 1:
            df = pd.DataFrame(self.df[keys[0]].copy())
        if len(keys) > 1:  # we can merge things
            for key in keys[1:]:
                df = pd.merge(df, pd.DataFrame(self.df[key]), on=["name", "url"])

        # For the quality_control pipeline
        columns = []
        for this in [
            "name",
            "url",
            "N_raw",
            "GC_raw_(%)",
            "Mean_quality_raw",
            "Phix_content_(%)",
            "Adapters_content_(%)",
            "Trimmed_reads_(%)",
            "N_final",
        ]:
            if this in df.columns:
                columns.append(this)
        df = df[columns]
        df.rename(columns={"name": "Sample name"}, inplace=True)

        from sequana.utils.datatables_js import DataTable

        datatable = DataTable(df, "multi_summary")
        datatable.datatable.datatable_options = {
            "scrollX": "300px",
            "pageLength": 30,
            "scrollCollapse": "true",
            "dom": "rtpB",
            "paging": "false",
            "buttons": ["copy", "csv"],
        }

        datatable.datatable.set_links_to_column("url", "Sample name")
        js = datatable.create_javascript_function()
        html_tab = datatable.create_datatable(float_format="%.3g")
        html = "{} {}".format(html_tab, js)

        self.jinja[
            "canvas"
        ] += """
    function onClick(e){
        window.open(e.dataPoint.url)
    }
}</script>"""

        caption = """<p>The table below gives a brief summary of the analysis. The
first column contains clickable sample name that redirects to complete summary
page. The table contains the following columns:</p>

   <b>Table caption</b>
    <table>
        <tr><td>N_raw</td><td>Number of reads in the data</td></tr>
        <tr><td>GC_raw_(%)</td><td>GC content in percentage in the raw data 
across all reads</td></tr>
        <tr><td>Mean_quality_raw</td><td>Mean quality across all reads all bases
in the raw data</td></tr>
        <tr><td>Phix_content_(%)</td><td>Percentage of reads found with Phix174</td></tr>
        <tr><td>Adapters_content_(%)</td><td>Percentage of reads with adapters (after phix
removal if applied)  </td></tr>
        <tr><td>Trimmed_reads_(%)</td><td>Percentage of reads trimmed (after
phix and adapter removal)</td></tr>
        <tr><td>N_final</td><td>Final number of reads (after phix and adapter
removal and trimming)</td></tr>
    </table>
"""
        infohtml = self.create_hide_section(
            "information", "(Show information)", caption, True
        )
        infohtml = "\n".join(infohtml)

        self.intro = introhtml + """ <hr><b>Summary</b>: """ + infohtml + html

        self.sections.append(
            {
                "name": None,
                "anchor": None,
                "content": self.jinja["canvas"] + "\n".join(self.jinja["sections"]),
            }
        )

    def _get_div(self, name, title):
        div = "<h2>%s</h2>" % title
        div += '<div id="chartContainer' + name
        div += '" style="height: 300px; width: 90%;"></div></hr>'
        return div

    def _get_df(self, method_name):
        data = getattr(self, method_name)()
        df = pd.DataFrame(
            {"name": self.get_projects(), "value": data, "url": self.get_urls()}
        )
        return df

    def populate_output_total_reads(self):
        df = self._get_df("get_output_total_reads")
        self.df["N_final"] = df.copy()
        self.df["N_final"].rename({"value": "N_final"}, axis=1, inplace=True)

    def populate_adapters(self):
        title = "Adapters content"
        df = self._get_df("get_adapters_percent")
        self.df["Adapters"] = df.copy()
        self.df["Adapters"].rename(
            {"value": "Adapters_content_(%)"}, axis=1, inplace=True
        )
        cb = CanvasBar(df, "Adapters content", "adapters", xlabel="Percentage")
        self.jinja["canvas"] += cb.to_html()
        self.jinja["sections"].append(self._get_div("adapters", title))

    def populate_nreads_raw(self):
        title = "Number of reads"
        df = self._get_df("get_nreads_raw")
        self.df["N_raw"] = df.copy()
        self.df["N_raw"].rename({"value": "N_raw"}, axis=1, inplace=True)
        cb = CanvasBar(
            df, "Number of reads (raw data)", "nreads_raw", xlabel="Number of reads"
        )
        self.jinja["canvas"] += cb.to_html()
        self.jinja["sections"].append(self._get_div("nreads_raw", title))

    def populate_mean_quality(self):
        title = "Mean quality (raw data)"
        df = self._get_df("get_mean_quality_samples")
        self.df["Mean_quality_raw"] = df.copy()
        self.df["Mean_quality_raw"].rename(
            {"value": "Mean_quality_raw"}, axis=1, inplace=True
        )
        cb = CanvasBar(df, title, "mean_quality", xlabel="mean quality")
        self.jinja["canvas"] += cb.to_html(options={"maxrange": 40})
        self.jinja["sections"].append(self._get_div("mean_quality", title))

    def populate_gc_samples(self):
        title = "GC content (raw)"
        df = self._get_df("get_gc_content_samples")
        self.df["GC_raw"] = df.copy()
        self.df["GC_raw"].rename({"value": "GC_raw_(%)"}, axis=1, inplace=True)
        cb = CanvasBar(df, title, "populate_gc_samples", xlabel="Percentage")
        self.jinja["canvas"] += cb.to_html(options={"maxrange": 100})
        self.jinja["sections"].append(self._get_div("populate_gc_samples", title))

    def populate_phix(self):
        title = "Phix content"
        df = self._get_df("get_phix_percent")
        self.df["Phix"] = df.copy()
        self.df["Phix"].rename({"value": "Phix_content_(%)"}, inplace=True, axis=1)
        cb = CanvasBar(df, title, "phix", xlabel="Percentage")
        self.jinja["canvas"] += cb.to_html()
        self.jinja["sections"].append(self._get_div("phix", title))

    def populate_trimming(self):
        title = "Trimming (raw data)"
        df = self._get_df("get_trimming_percent")
        self.df["Trimmed"] = df.copy()
        self.df["Trimmed"].rename({"value": "Trimmed_reads_(%)"}, axis=1, inplace=True)
        cb = CanvasBar(df, title, "trimming", xlabel="Percentage")
        self.jinja["canvas"] += cb.to_html()
        self.jinja["sections"].append(self._get_div("trimming", title))

    def get_projects(self):
        return [x.data["project"] for x in self.summaries]

    def get_urls(self):
        return [
            x.replace("summary.json", "summary.html") for x in self.get_unique_names()
        ]

    def get_unique_names(self):
        """reduce the filenames length removing the common suffix
        remove also the filename itself.
        """
        projects = self.get_projects()
        if len(projects) == set(projects):
            return projects
        else:
            return self.filenames

    ##########################################################################
    # The retrieval of all specific data for all summaries

    def get_adapters_percent(self):
        return [x.get_adapters_percent() for x in self.summaries]

    def get_nreads_raw(self):
        return [x.get_nreads_raw() for x in self.summaries]

    def get_phix_percent(self):
        return [x.get_phix_percent() for x in self.summaries]

    def get_gc_content_samples(self):
        return [x.get_fastq_stats_samples()["GC content"]["R1"] for x in self.summaries]

    def get_trimming_percent(self):
        return [x.get_trimming_percent() for x in self.summaries]

    def get_mean_quality_samples(self):
        return [x.get_mean_quality_samples() for x in self.summaries]

    def get_output_total_reads(self):
        return [x.get_output_total_reads() for x in self.summaries]

    def parse(self):
        pass
