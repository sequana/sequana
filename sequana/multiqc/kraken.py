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
log = logging.getLogger('multiqc.sequana/kraken')


class MultiqcModule(BaseMultiqcModule):
    """

    Here keys are the sample name

    {   
     "25": {
      "Bacteria": 0.4,
      "Metazoa": 2.0,
      "Unclassified": 0.0,
      "Viruses": 97.6
     },
     "26": {
      "Bacteria": 8.800000000000002,
      "Metazoa": 0.4,
      "Unclassified": 18.4,
      "Viruses": 72.39999999999998
     },
    }

    .. warning:: this is not really a multiqc here because we read a file
        that summarizes several samples together. 

    """
    def __init__(self):

        # Initialise the parent object
        super(MultiqcModule, self).__init__(
            name='Sequana/Kraken',    # name that appears at the top
            anchor='sequana_kraken',  #
            target='sequana',  # Name show that link to the following href
            href='http://github.com/sequana/sequana/',
            info="sequana_kraken pipeline multi summary")

        self.sequana_data = {}

        # In theory only one file
        for myfile in self.find_log_files("sequana_kraken"):
            name = myfile['s_name']
            d = json.loads(myfile['f'])


            for k,v in d.items():
                # check for None and NAN
                S = sum(list(v.values()))
                U = v['Unclassified']
                v['Classified'] = S - U
                self.sequana_data[k] = v

        if len(self.sequana_data) == 0:
            log.debug("No samples found: sequana_kraken")
            raise UserWarning


        """info = "<ul>"
        for this in sorted(self.sequana_data.keys()):
            info += '<li><a href="{}/summary.html">{}</a></li>'.format(this,this,this)
        info += "</ul>"
        href="http://sequana.readthedocs.io/en/master/"
        target = "Sequana"
        mname = '<a href="{}" target="_blank">{}</a> individual report pages:'.format(href, target)
        self.intro = '<p>{} {}</p>'.format( mname, info)
        """

        log.info("Found {} reports".format(len(self.sequana_data)))

        self.populate_columns()
        self.add_kraken()

    def add_classified_vs_unclassified(self):
        pass

    def _set_nan_to_zero(self, x):
        try:
            x + 0
            return x
        except:
            return 0

    def add_kraken(self):
        data = {}
        

        for name in self.sequana_data.keys():

            for this in ['Viruses', 'Bacteria', 'Unclassified', 'Metazoa']:
                if this not in self.sequana_data[name]:
                    self.sequana_data[name][this] = 0
            data[name] = {
                'viruses': self._set_nan_to_zero(self.sequana_data[name]['Viruses']),
                'bacteria': self._set_nan_to_zero(self.sequana_data[name]['Bacteria']),
                'unclassified': self._set_nan_to_zero(self.sequana_data[name]['Unclassified']),
                'metazoa': self._set_nan_to_zero(self.sequana_data[name]['Metazoa'])
            }

        pconfig = {
            "title": "Taxonomy",
            "percentages": True,
            "min": 0,
            "max":100,
            "format": '{0:.2f}',
            "logswitch": False,
        }

        keys = OrderedDict()
        keys['viruses'] = {'color': '#437bb1', 'name': 'Viruses'}
        keys['bacteria'] = {'color': '#b1084c', 'name': 'Bacteria'}
        keys['metazoa'] = {'color': 'green', 'name': 'Metazoa'}
        keys['unclassified'] = {'color': 'grey', 'name': 'Unclassified'}

        self.add_section(
            name = 'taxonomy',
            anchor = 'taxonomy',
            description = 'The following barplots summarizes the kraken analysis for each sample. ',
            helptext = "",
            plot = bargraph.plot(data, keys, pconfig))

    def populate_columns(self):
        headers = {}

        if any(['Classified' in self.sequana_data[s] for s in self.sequana_data]):
            headers['Classified'] = {
                'title': 'Classified reads (%)',
                'description': 'classified reads (%)',
                'min': 0,
                'max': 100,
                'scale': 'RdYlGn',
                'format': '{0:.2f}',
                'shared_key': 'count',
            }

        if len(headers.keys()):
            self.general_stats_addcols(self.sequana_data, headers)

