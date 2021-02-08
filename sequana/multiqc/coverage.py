#!/usr/bin/env python

""" MultiQC module to parse output from sequana (coverage)"""
import os
import re
from math import log10

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

    def __init__(self):

        # Initialise the parent object
        super(MultiqcModule, self).__init__(
            name='Sequana/coverage',    # name that appears at the top
            anchor='sequana_coverage',
            target='sequana',  # Name show that link to the following href
            href='http://github.com/sequana/sequana/',
            info="sequana_coverage multi Summary")

        self.sequana_data = dict()
        self.sequana_desc = {}
        self.sequana_root = {}
        self.sequana_caller = {}

        # Do we have a unique chromosome/contig name ?
        chrom_names = set([])
        for myfile in self.find_log_files("sequana_coverage"):
            name = myfile['s_name']
            if name.startswith("summary_"):
                name = name.replace("summary_", "")
            data = self.parse_logs(myfile["f"])
            chrom_names.add(data['data']['chrom_name'])

        self._maps = {}
        if len(chrom_names) == 1:
            self.unique = True
        else:
            self.unique = False

        for myfile in self.find_log_files("sequana_coverage"):
            name = myfile['s_name']
            if name.startswith("summary_"):
                name = name.replace("summary_", "")
            data = self.parse_logs(myfile["f"])
            if self.unique is True:
                 key = data["sample_name"]
            else:
                # convert into string in case
                key = data["sample_name"]+"/"+ str(data['data']['chrom_name'])

            # in any case, just keep track of key / chrom_name

            # sometimes keys/name of chromosomes are just integer, so we need to
            # convert to string
            self._maps[key] = data["sample_name"] + "/" + str(data['data']['chrom_name'])

            self.sequana_data[key] = data['data']
            self.sequana_desc[key] = data['data_description']
            self.sequana_root[key] = myfile['root']
            try:
                self.sequana_caller[key] = data['caller']
            except:
                self.sequana_caller[key] = "undefined"
        if len(self.sequana_data) == 0:
            log.debug("No samples found: sequana_coverage")
            raise UserWarning

        info = "<ul>"
        for this in sorted(self.sequana_data.keys()):
            sample_name, chrom_name = self._maps[this].split("/")
            #sample_name, chrom_name = this.split("/")

            root = self.sequana_root[this]

            # sequana_coverage store multiqc locally
            # the pipelines stores it in ./multiqc/ the only way to have the
            # correct path is to figure out who creates the summary file for
            # multiqc. If sequana_coverage, no need to alter the path but if
            # the multiqc were create by the pipelines, you need to come back to
            # the root of the pipeline
            if self.sequana_caller[this] in ["sequana_coverage",
                    "sequana.bedtools"]:
                subpath = ""
            else:
                subpath = "../"

            info += '<li><a href="{}/{}.cov.html">{}</a></li>'.format(
                subpath + root,
                chrom_name,  this)
        info += "</ul>"
        href="http://sequana.readthedocs.io/en/master/"
        target = "Sequana"
        mname = '<a href="{}" target="_blank">{}</a> individual report pages:'.format(href, target)
        self.intro = '<p>{} {}</p>'.format( mname, info)

        log.info("Found {} reports".format(len(self.sequana_data)))

        self.populate_columns()
        self.add_hist_coverage()
        self.add_DOC()
        self.add_BOC()
        try:self.add_CV()
        except:pass
        self.add_length()
        self.add_ROI()
        self.add_C3()

    def parse_logs(self, log_dict):
        import json
        log_dict = json.loads(log_dict)
        return log_dict

    def add_BOC(self):
        data = {}
        for name in self.sequana_data.keys():
            data[name] = {'BOC': self.sequana_data[name]["BOC"]}

        pconfig = {
            "title": self.sequana_desc[name]["BOC"],
            "percentages": False,
            "min": 0,
            "max": 0,
            "logswitch": False}
        self.add_section(
            name = 'Breadth of coverage',
            anchor = 'boc',
            description = 'Breadth of coverage: proportion of the genome ' + \
                'covered by at least one read',
            helptext = "",
            plot = bargraph.plot(data, None, pconfig))

    def add_DOC(self):
        data = {}
        for name in self.sequana_data.keys():
            data[name] = {'DOC': self.sequana_data[name]["DOC"]}

        pconfig = {
            "title": self.sequana_desc[name]["DOC"],
            "percentages": False,
            "min": 100,
            "logswitch": True}
        self.add_section(
            name = 'Depth of coverage',
            anchor = 'doc',
            description = 'Depth of coverage: average number of reads' + \
                ' mapping on each genome position',
            helptext = "",
            plot = bargraph.plot(data, None, pconfig))

    def add_CV(self):
        data = {}
        for name in self.sequana_data.keys():
            data[name] = {'CV': self.sequana_data[name]["CV"]}

        pconfig = {
            "title": self.sequana_desc[name]["CV"],
            "percentages": False,
            "min": 100,
            "logswitch": False}

        self.add_section(
            name = 'Coefficient of Variation',
            anchor = 'cv',
            description = 'The ratio of DOC mean by DOC standard deviation',
            helptext = "",
            plot = bargraph.plot(data, None, pconfig))

    def add_ROI(self):
        data = {}
        for name in self.sequana_data.keys():
            data[name] = {'ROI': self.sequana_data[name]["ROI"]}

        pconfig = {
            "title": self.sequana_desc[name]["ROI"],
            "percentages": False,
            "min": 100,
            "logswitch": False}

        self.add_section(
            name = 'ROI',
            anchor = 'ROI',
            description = 'Number of regions of interest',
            helptext = "",
            plot = bargraph.plot(data, None, pconfig))

    def add_C3(self):
        data = {}
        for name in self.sequana_data.keys():
            data[name] = {'C3': self.sequana_data[name]["C3"]}

        pconfig = {
            "title": self.sequana_desc[name]["C3"],
            "percentages": False,
            "min": 0,
            "max": 0,
            "logswitch": False}

        self.add_section(
            name = 'C3',
            anchor = 'C3',
            description = 'Centralness (roughly speaking, ratio of ' + \
                          'outliers versus total genome length).',
            helptext = "",
            plot = bargraph.plot(data, None, pconfig))

    def add_length(self):
        data = {}
        for name in self.sequana_data.keys():
            data[name] = {'length': self.sequana_data[name]["length"]}

        pconfig = {
            "title": self.sequana_desc[name]["length"],
            "percentages": False,
            "min": 0,
            "logswitch": False}
        self.add_section(
            name = 'Contig length',
            anchor = 'length',
            description = 'Length of the contig/chromosome',
            helptext = "",
            plot = bargraph.plot(data, None, pconfig))

    def add_hist_coverage(self):
        data = dict()
        data_norm = dict()
        for s_name in self.sequana_data:
            try:
                X = self.sequana_data[s_name]["hist_coverage"]['X']
                Y = self.sequana_data[s_name]["hist_coverage"]['Y']
                Y = [y  for y in Y]
                data[s_name] = {x: y if y else 0 for x,y in zip(X, Y)}
            except KeyError:
                pass

        if len(data) == 0:
            log.debug('no data for the coverage plots')
            return None

        pconfig = {
            'id': 'sequana_coverage_hist',
            'title': 'Depth of Coverage',
            'ylab': '#',
            'xlab': 'DOC',
            'ymin': 0,
            #'ymax': 0.2,
            #'xmax': 100,
            'xmin': 0,
            'yDecimals': True,
            'tt_label': '<b>{point.x:.2f} DOC</b>: {point.y:.4f}',
            #'colors': self.get_status_cols('per_sequence_gc_content'),
            'data_labels': [
                #{'name': 'Percentages', 'ylab': 'Percentage'},
                {'name': 'Counts', 'ylab': 'PDF'}
            ]
        }

        self.add_section (
            name = 'Depth of Coverage Histogram',
            anchor = 'coverage_hist',
            description = ("Histogram (normalised) of the depth of coverage. For"
                        " convenience, only  99% the data (centered) to "
                        "avoid outliers. For detailled histograms, please see "
                        " the links above "),
            #plot = linegraph.plot([data_norm, data], pconfig))
            plot = linegraph.plot(data, pconfig))


    def populate_columns(self):
        headers = {}
        formats = {
            "BOC": '{:,.2f}',
            "CV": '{:,.2f}',
            "DOC": '{:,.2f}',
            "length": '{:,d}',
            "ROI": '{:,d}',
            "C3": '{:,.2f}'
        }

        for field in ['BOC', 'DOC', 'ROI', 'length', "CV", "C3"]:

            # description are supposed to be the same for all samples, let us
            # take the first one
            desc = [self.sequana_desc[s][field] for s in self.sequana_desc][0]

            if any([field in self.sequana_data[s] for s in self.sequana_data]):
                headers[field] = {
                    'title': field,
                    'description': desc,
                    'min':0,
                    'max':0,
                    'scale': 'RdYlGn',
                    'format': formats[field],
                    'shared_key': 'count',
                }

                if field in ['BOC']:
                    headers[field]['min'] = 0
                    headers[field]['max'] = 100


        if len(headers.keys()):
            self.general_stats_addcols(self.sequana_data, headers)


