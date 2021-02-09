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
import colorlog
logger = colorlog.getLogger(__name__)



class MultiqcModule(BaseMultiqcModule):
    """

    """
    def __init__(self):

        # Initialise the parent object
        super(MultiqcModule, self).__init__(
            name='Sequana/pacbio_amplicon',    # name that appears at the top
            anchor='sequana_pacbio_amplicon',  #
            target='sequana',  # Name show that link to the following href
            href='http://github.com/sequana/sequana/',
            info="amplicon pipeline multi summary")

        self.sequana_data = {}

        # In theory only one file
        for myfile in self.find_log_files("sequana_pacbio_amplicon"):
            logging.info("Parsing {}".format(myfile))
            name = myfile['s_name']
            name = name.replace("sequana_pacbio_amplicon_", "")
            self.sequana_data[name] = self.parse_logs(myfile["f"])

        if len(self.sequana_data) == 0:
            logger.debug("No samples found: sequana_pacbio_amplicon")
            raise UserWarning


        logger.info("Found {} reports".format(len(self.sequana_data)))

        self.populate_columns()
        self.add_ccs_reads()
        self.add_ccs_reads_hist()

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
            name = 'CCS reads',
            anchor = 'ccs_reads',
            description = 'The following barplots summarizes the number of CCS reads.',
            helptext = "",
            plot = bargraph.plot(data, None, pconfig))

    def add_ccs_reads_hist(self):
        data = {}
        #for name in self.sequana_data:
        #    data = [item['ccs_reads'] for item in self.sequana_data.items()]

        
        data = [self.sequana_data[key]["ccs_reads"] for key in self.sequana_data.keys()]
        from pylab import hist
        Y, X, _ = hist(data, bins=20)

        data = {}
        data['sample'] = {x:y for x,y in zip(X, Y)}

        pconfig = {
            "title": "CCS reads ",
            "percentages": False,
            "min": 0,
            #"max":100,
            "format": '{0:.2f}'
        }

        self.add_section(
            name = 'CCS reads histogram',
            anchor = 'ccs_reads_hist',
            description = 'CCS reads histogram.',
            helptext = "",
            plot = linegraph.plot(data, pconfig))

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

