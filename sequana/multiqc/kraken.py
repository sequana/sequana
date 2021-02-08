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

import colorlog
log = colorlog.getLogger(__name__)



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
        self.add_classified_vs_unclassified()

    def add_classified_vs_unclassified(self):
        data = {}
        for sample_name in self.sequana_data.keys():
            C = self.sequana_data[sample_name]['Classified']
            U = self.sequana_data[sample_name]['Unclassified']
            data[sample_name] = {'unclassified': U, 'classified': C}

        pconfig = {
            "title": "classification",
            "cpswitch": False,
            "min": 0,
            "max":100,
            "format": '{0:.2f}',
            "logswitch": False,
        }

        keys = OrderedDict()
        keys['classified'] = {'color': 'green', 'name': 'Classified'}
        keys['unclassified'] = {'color': 'red', 'name': 'Unclassified'}

        self.add_section(
            name = 'Classification ratio',
            anchor = 'classification_ratio',
            description = 'The following barplots shows proportion of classified and unclassified reads in percentage',
            helptext = "",
            plot = bargraph.plot(data, keys, pconfig))

    def _set_nan_to_zero(self, x):
        try:
            x + 0
            return x
        except:
            return 0

    def add_kraken(self):
        data = {}

        # First, we figure out all possible names
        kingdoms = set([x for k in self.sequana_data.keys() for x in self.sequana_data[k].keys()])

        colors = ['Archaea', 'Bacteria', 'Eukaryota', 'Viruses', 'Metazoa',
            'Fungi', "Unclassified", "Classified"]

        for sample_name in self.sequana_data.keys():
            for kingdom in sorted(kingdoms):
                if kingdom not in self.sequana_data[sample_name]:
                    self.sequana_data[sample_name][kingdom] = 0

            data[sample_name] = {"others": 0}
            for kingdom in sorted(kingdoms):

                if kingdom not in colors:
                    # here we add together non-superkingdom + other artifical
                    # sequences
                    data[sample_name]["others"] += \
                         self._set_nan_to_zero(self.sequana_data[sample_name][kingdom])
                else:
                    data[sample_name][kingdom.lower()] = \
                         self._set_nan_to_zero(self.sequana_data[sample_name][kingdom])
            data[sample_name]['unclassified'] = \
                self._set_nan_to_zero(self.sequana_data[sample_name]['Unclassified'])

        pconfig = {
            "title": "Taxonomy by kingdom",
            #"percentages": True,
            "cpswitch": False,
            "min": 0,
            "max":100,
            "format": '{0:.2f}',
            "logswitch": False,
        }

        keys = OrderedDict()
        # superkingdom:
        keys['archea'] = {'color': 'orange', 'name': 'Archea'}
        keys['bacteria'] = {'color': '#b1084c', 'name': 'Bacteria'}
        keys['eukaryota'] = {'color': 'green', 'name': 'Eukaryota'}
        keys['viruses'] = {'color': '#437bb1', 'name': 'Viruses'}
        # kingdom: 
        keys['metazoa'] = {'color': 'green', 'name': 'Metazoa'}
        keys['fungi'] = {'color': 'purple', 'name': 'Fungi'}
        # others
        keys['unclassified'] = {'color': 'grey', 'name': 'Unclassified'}
        keys['others'] = {'color': 'blue', 'name': 'Others'}
        # subkingdom
        #keys['viridiplantae'] = {'color': 'yellow', 'name': 'Viridiplantae'}
        #keys['dikarya'] = {'color': 'brown', 'name': 'dikarya'}

        self.add_section(
            name = 'Taxonomy by kingdom',
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

        for name in ['Viruses', 'Bacteria', 'Eukaryota', 'Archea', 'Fungi']:
            if any([name in self.sequana_data[s] for s in self.sequana_data]):
                headers[name] = {
                    'title': 'Reads classified as {} (%)'.format(name.lower()),
                    'description': 'Reads classified as {} (%)'.format(name.lower()),
                    'min': 0,
                    'max': 100,
                    'scale': 'RdYlGn',
                    'format': '{0:.2f}',
                    'shared_key': 'count',
                }

        if len(headers.keys()):
            self.general_stats_addcols(self.sequana_data, headers)

