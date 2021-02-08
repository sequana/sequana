#!/usr/bin/env python
""" MultiQC module to parse output from sequana"""
import os
import re

import logging
import json

logging.captureWarnings(True)
from multiqc import config
from multiqc.plots import linegraph, table, heatmap, bargraph
from multiqc.modules.base_module import BaseMultiqcModule
logging.captureWarnings(False)

# Initialise the logger
import colorlog
log = colorlog.getLogger(__name__)



class MultiqcModule(BaseMultiqcModule):

    def __init__(self):

        # Initialise the parent object
        super(MultiqcModule, self).__init__(
            name='Sequana/quality_control',
            anchor='sequana_quality_control',
            target='sequana_quality_control',
            href='http://github.com/sequana/sequana/',
            info="(sequana quality control multi summary)")

        self.data = {}
        for myfile in self.find_log_files("sequana_quality_control"):
            thisdata =  self.parse_logs(myfile["f"])
            name = thisdata["project"]
            self.data[name] = self.parse_logs(myfile["f"])

        if len(self.data) == 0:
            log.debug("No samples found: sequana_quality_control")
            raise UserWarning

        log.info("Found {} reports".format(len(self.data)))

        name = list(self.data.keys())[0]
        if "phix_section" in self.data[name]:
            self._has_phix = True
        else:
            self._has_phix = False


        self.populate_columns()
        self.add_phix_section()
        self.add_gc_section()
        self.add_adapter_section()

    def populate_columns(self):

        headers = {}
        if any(['multiqc_total_reads' in self.data[s] for s in self.data]):
            headers['multiqc_total_reads'] = {
                'title': 'number of reads',
                'description': 'number of reads',
                #'max': 100,
                'min': 0,
                #'modify': lambda x: x * 100,
                'scale': 'RdYlGn',
                #'format': '{:,.1f}'
                'shared_key': 'multiqc_total_reads',
                'hidden': True,
            }
        if self._has_phix:
            headers['phix'] = {
                    'title': 'Phix', 
                    'min':0,
                    'max':100,
                    'hidden': True
                }

        if len(headers.keys()):
            self.general_stats_addcols(self.data, headers)

    def parse_logs(self, log_dict):
        data = json.loads(log_dict)
        this =  data["cutadapt_json"]["Number of reads"]["Total paired reads"]
        data["multiqc_total_reads"] = this
        # Phix removal is optional
        if "phix_section" in data.keys():
            this =  data["phix_section"]["contamination"]
            data["phix"] = this
        return data

    def add_gc_section(self):
        data = {}
        for name in self.data.keys():
            data[name] = {'GC': self.data[name]["fastq_stats_samples_json"]["GC content"]['R1']}
        pconfig = {
            "title": "Percentage of GC in the raw dat",
            "percentages": True,
            "min": 0,
            "max": 100
        }

        self.add_section(
            name = 'GC percent',
            anchor = 'GC percent',
            description = "",
            helptext = "",
            plot = bargraph.plot(data, None, pconfig)
        )

    def add_phix_section(self):
        pconfig = {
            "title": "Percentage of phix in the raw data",
            "percentages": True,
            "min": 0,
            "max": 100
        }
        if self._has_phix:
            data = {}
            for name in self.data.keys():
                data[name] = {'phix_qc': self.data[name]["phix_section"]["contamination"]}
            description = 'Amount of Phix present in the raw data (and removed)'
            if sum([float(x['phix_qc']) for x in data.values()]) == 0:
                plot = None
                description += "\n\nNo Phix found in the data."
            else:
                plot = bargraph.plot(data, None, pconfig)
        else:
            description = "No Phix removed (there may be some)"
            plot = None


        self.add_section(
            name = 'Phix presence',
            anchor = 'mean_read_length',
            description = description,
            helptext = "",
            plot = plot)

    def add_adapter_section(self):
        data = {}
        for name in self.data.keys():
            thisdata = self.data[name]["cutadapt_json"]["percent"]["Pairs kept"]
            data[name] = {'pairs_kept': thisdata.replace("(","").replace(")","").replace("%","")}

        pconfig = {
            "title": "Percentage of pairs kept",
            "percentages": True,
            "min": 0,
            "max": 100,
        }

        self.add_section(
            name = 'Pairs kept',
            anchor = 'Pairs kept',
            description = 'Pairs kept',
            helptext = "",
            plot = bargraph.plot(data, None, pconfig))







