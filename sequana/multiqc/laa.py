#!/usr/bin/env python

""" MultiQC module to parse output from sequana"""
import os
import re
import json

from collections import OrderedDict
# prevent boring warning (version 1.0)
import logging
logging.captureWarnings(True)
from multiqc import config
from multiqc.modules.base_module import BaseMultiqcModule
from multiqc.plots import linegraph, table, heatmap, bargraph
logging.captureWarnings(False)

# Initialise the logger
log = logging.getLogger('multiqc.sequana/laa')


class MultiqcModule(BaseMultiqcModule):
    """

    """
    def __init__(self):

        # Initialise the parent object
        super(MultiqcModule, self).__init__(
            name='Sequana/laa',    # name that appears at the top
            anchor='sequana_laa',  #
            target='sequana',  # Name show that link to the following href
            href='http://github.com/sequana/sequana/',
            info="sequana laa pipeline multi summary")

        self.sequana_data = {}

        # In theory only one file
        for myfile in self.find_log_files("sequana_laa"):
            logging.info("Parsing {}".format(myfile))
            name = myfile['s_name']
            name = name.replace("sequana_laa_", "")
            self.sequana_data[name] = self.parse_logs(myfile["f"])

        print(self.sequana_data)

        if len(self.sequana_data) == 0:
            log.debug("No samples found: sequana_laa")
            raise UserWarning


        log.info("Found {} reports".format(len(self.sequana_data)))

        self.populate_columns()
        self.add_ccs_reads()


    def parse_logs(self, log_dict):
        import json
        log_dict = json.loads(log_dict)
        return log_dict

    def add_ccs_reads(self):
        data = {}

        for name in self.sequana_data.keys():
            data[name] = {
                'ccs_reads': self.sequana_data[name]['ccs_reads']
            }

        pconfig = {
            "title": "CCS reads",
            "percentages": False,
            "min": 0,
            #"max":100,
            "format": '{0:.2f}',
            "logswitch": False,
        }

        self.add_section(
            name = 'taxonomy',
            anchor = 'taxonomy',
            description = 'The following barplots summarizes the kraken analysis for each sample. ',
            helptext = "",
            plot = bargraph.plot(data, None, pconfig))

    def populate_columns(self):
        headers = {}

        if any(['ccs_reads' in self.sequana_data[s] for s in self.sequana_data]):
            headers['ccs_reads'] = {
                'title': 'CCS reads ',
                'description': 'CCS reads',
                'min': 0,
                #'max': 100,
                'scale': 'RdYlGn',
                'format': '{0:.2f}',
                'shared_key': 'count',
            }

        if len(headers.keys()):
            self.general_stats_addcols(self.sequana_data, headers)

