#!/usr/bin/env python

""" MultiQC module to parse output from sequana"""
import os
import re

# prevent boring warning (version 1.0)
import logging
logging.captureWarnings(True)
from multiqc import config
from multiqc.modules.base_module import BaseMultiqcModule
from multiqc.plots import linegraph, table, heatmap, bargraph
logging.captureWarnings(False)

# Initialise the logger
log = logging.getLogger('multiqc.sequana/bamtools_stats')


class MultiqcModule(BaseMultiqcModule):

    def __init__(self):

        # Initialise the parent object
        super(MultiqcModule, self).__init__(
            name='Sequana/bamtools_stats',    # name that appears at the top
            anchor='sequana',  #
            target='sequana',  # Name show that link to the following href
            href='http://github.com/sequana/sequana/',
            info="pipelines multi Summary")

        self.sequana_data = {}
        for myfile in self.find_log_files("sequana_bamtools_stats"):
            logging.info("Parsing {}".format(myfile))
            #print( myfile['f'] )       # File contents
            #print( myfile['s_name'] )  # Sample name (from cleaned filename)
            #print( myfile['fn'] )      # Filename
            #print( myfile['root'] )    # Directory file was in
            name = myfile['s_name']
            if name.startswith("sequana_bamtools_stats_"):
                name = name.replace("sequana_bamtools_stats_", "")
            self.sequana_data[name] = self.parse_logs(myfile["f"])

        """info = "<ul>"
        for this in sorted(self.sequana_data.keys()):
            info += '<li><a href="{}/summary.html">{}</a></li>'.format(this,this,this)
        info += "</ul>"
        href="http://sequana.readthedocs.io/en/master/"
        target = "Sequana"
        mname = '<a href="{}" target="_blank">{}</a> individual report pages:'.format(href, target)
        self.intro = '<p>{} {}</p>'.format( mname, info)
        """

        if len(self.sequana_data) == 0:
            log.debug("Could not find any data in {}".format(config.analysis_dir))
            raise UserWarning

        log.info("Found {} reports".format(len(self.sequana_data)))

        self.populate_columns()
        self.add_total_read_section()
        self.add_mapped_reads_section()
        self.add_duplicates_section()
        self.add_forward_reverse_section()
        self.add_failed_QC()

    def parse_logs(self, log_dict):
        """Parse this kind of logs::

            **********************************************
            Stats for BAM file(s):
            **********************************************

            Total reads:       496
            Mapped reads:      491  (98.9919%)
            Forward strand:    235  (47.379%)
            Reverse strand:    261  (52.621%)
            Failed QC:         0    (0%)
            Duplicates:        0    (0%)
            Paired-end reads:  0    (0%)
        """
        data = {}
        for line in log_dict.split("\n"):
            if line.strip() == "" or line.startswith("*") or "Stats for BAM" in line:
                continue
            key, value = line.split(":", 1)
            data[key] = float(value.split()[0])
        return data

    def add_failed_QC(self):
        data = {}
        for name in self.sequana_data.keys():
            data[name] = {'failed_qc': self.sequana_data[name]["Failed QC"]}

        pconfig = {
            "title": "Failed QC (%)",
            "percentages": True,
            "min": 0,
            "max":100,
            "logswitch": False,
        }

        self.add_section(
            name = 'Failed QC',
            anchor = 'failed_qc',
            description = 'Failed QC',
            helptext = "",
            plot = bargraph.plot(data, None, pconfig))

    def add_forward_reverse_section(self):
        data = {}
        for name in self.sequana_data.keys():
            forward = self.sequana_data[name]["Forward strand"]
            reverse = self.sequana_data[name]["Reverse strand"]
            data[name] = {'fwd_rev': forward}

        pconfig = {
            "title": "Forward/reverse",
            "percentages": True,
            "min": 0,
            "max":100,
            "logswitch": False,
        }

        self.add_section(
            name = 'Forward/Reverse',
            anchor = 'fwd_rev',
            description = 'Forward reverse',
            helptext = "",
            plot = bargraph.plot(data, None, pconfig))


    def add_total_read_section(self):
        data = {}
        for name in self.sequana_data.keys():
            data[name] = {'total_read': self.sequana_data[name]["Total reads"]}

        pconfig = {
            "title": "Total reads",
            "percentages": False,
            "min": 0,
            "logswitch": True,
        }

        self.add_section(
            name = 'Total reads',
            anchor = 'total_read',
            description = 'Total reads found in the BAM file.',
            helptext = "",
            plot = bargraph.plot(data, None, pconfig))

    def add_mapped_reads_section(self):
        data = {}
        for name in self.sequana_data.keys():
            total = self.sequana_data[name]["Total reads"]
            mapped_reads = self.sequana_data[name]["Mapped reads"]
            data[name] = {'mapped_reads': 100 * mapped_reads / total}

        pconfig = {
            "title": "Mapped reads",
            "percentages": False,
            "min": 0,
            "max": 100,
            "logswitch": True,
        }

        self.add_section(
            name = 'Mapped reads',
            anchor = 'mapped reads',
            description = 'Total mapped reads found in the BAM file.',
            helptext = "",
            plot = bargraph.plot(data, None, pconfig))

    def add_duplicates_section(self):
        data = {}
        for name in self.sequana_data.keys():
            total = float(self.sequana_data[name]["Total reads"])
            duplicates = float(self.sequana_data[name]["Duplicates"])
            print(duplicates)
            data[name] = {'duplicates': 100 * duplicates/total}

        pconfig = {
            "title": "Duplicates",
            "percentages": True,
            "min": 0,
            "max": 100,
            "logswitch": False,
        }

        self.add_section(
            name = 'Duplicates',
            anchor = 'duplicates',
            description = 'Duplicated reads.',
            helptext = "",
            plot = bargraph.plot(data, None, pconfig))

    def populate_columns(self):
        headers = {}

        if any(['Total reads' in self.sequana_data[s] for s in self.sequana_data]):
            headers['Total reads'] = {
                'title': 'Read count',
                'description': 'read count',
                'min': 0,
                'scale': 'RdYlGn',
                'format': '{:,0d}',
                'shared_key': 'count',
            }


        if len(headers.keys()):
            self.general_stats_addcols(self.sequana_data, headers)

