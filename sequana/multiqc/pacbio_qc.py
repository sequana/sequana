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
import colorlog
logger = colorlog.getLogger(__name__)



class MultiqcModule(BaseMultiqcModule):

    def __init__(self):

        # Initialise the parent object
        super(MultiqcModule, self).__init__(
            name='Sequana/pacbio',    # name that appears at the top
            anchor='sequana',  # ??
            target='sequana',  # Name show that link to the following href
            href='http://github.com/sequana/sequana/',
            info="pipelines multi Summary")

        self.sequana_data = {}
        for myfile in self.find_log_files("sequana_pacbio_qc"):
            logging.info("Parsing {}".format(myfile))
            #print( myfile['f'] )       # File contents
            #print( myfile['s_name'] )  # Sample name (from cleaned filename)
            #print( myfile['fn'] )      # Filename
            #print( myfile['root'] )    # Directory file was in
            name = myfile['s_name']
            if name.startswith("sequana_summary_pacbio_qc_"):
                name = name.replace("sequana_summary_pacbio_qc_", "")
            self.sequana_data[name] = self.parse_logs(myfile["f"])

        if len(self.sequana_data) == 0:
            logger.debug("No samples found: sequana_pacbio_qc")
            raise UserWarning

        info = "<ul>"
        for this in sorted(self.sequana_data.keys()):
            info += '<li><a href="../{}/summary.html">{}</a></li>'.format(this,this)
        info += "</ul>"
        href="http://sequana.readthedocs.io/en/master/"
        target = "Sequana"
        mname = '<a href="{}" target="_blank">{}</a> individual report pages:'.format(href, target)
        self.intro = '<p>{} {}</p>'.format( mname, info)


        logger.info("Found {} reports".format(len(self.sequana_data)))

        self.populate_columns()
        self.add_count_section()
        self.add_read_length_section()
        self.add_hist_GC()
        self.add_hist_length()
        # Add histogram GC, SNR, read length

    def add_read_length_section(self):
        data = {}
        for name in self.sequana_data.keys():
            data[name] = {'mean': self.sequana_data[name]["mean_length"]}

        pconfig = {
            "title": "Mean read length",
            "percentages": False,
            "min": 100,
            "logswitch": True,
        }

        self.add_section(
            name = 'Mean read length',
            anchor = 'mean_read_length',
            description = 'Mean length of the reads per sample',
            helptext = "",
            plot = bargraph.plot(data, None, pconfig))

    def add_count_section(self):
        data = {}
        for name in self.sequana_data.keys():
            data[name] = {'count': self.sequana_data[name]["count"]}

        pconfig = {
            "title": "Number of reads per sample",
            "plot_tt_percentages": False,
            "min": 100,
            "logswitch": True,
        }

        self.add_section(
            name = 'Number of reads',
            anchor = 'number_of_reads',
            description = 'Number of reads per sample.',
            helptext = "",
            plot = bargraph.plot(data, None, pconfig)
            )

    def add_hist_GC(self):
        """ Create the HTML for the FastQC GC content plot """
        data = dict()
        data_norm = dict()
        for s_name in self.sequana_data:
            try:
                X = self.sequana_data[s_name]["hist_gc"]['X']
                Y = self.sequana_data[s_name]["hist_gc"]['Y']
                Y = [y / float(sum(Y)) for y in Y]
                data[s_name] = {x:10*y for x,y in zip(X[1:], Y)}
            except KeyError:
                pass
            #else:
            #    data_norm[s_name] = dict()
            #    total = sum( [ c for c in data[s_name].values() ] )
            #    for gc, count in data[s_name].items():
            #        data_norm[s_name][gc] = (count / total) * 100

        if len(data) == 0:
            logger.debug('no data for the GC content plots')
            return None

        pconfig = {
            'id': 'sequana_pacbio_per_sequence_gc_content_plot',
            'title': 'Per Sequence GC Content',
            'ylab': 'Count',
            'xlab': '% GC',
            'ymin': 0,
            #'ymax': 0.2,
            'xmax': 100,
            'xmin': 0,
            'yDecimals': False,
            'tt_label': '<b>{point.x}% GC</b>: {point.y}',
            #'colors': self.get_status_cols('per_sequence_gc_content'),
            'data_labels': [
                {'name': 'Percentages', 'ylab': 'Percentage'},
                {'name': 'Counts', 'ylab': 'PDF'}
            ]
        }

        self.add_section (
            name = 'Per Sequence GC Content',
            anchor = 'fastqc_per_sequence_gc_content',
            description = "GC content (normalised)",
            #plot = linegraph.plot([data_norm, data], pconfig))
            plot = linegraph.plot(data, pconfig))

    def add_hist_length(self):
        """ Create the HTML for the FastQC GC content plot """
        data = dict()
        data_norm = dict()
        for s_name in self.sequana_data:
            X = self.sequana_data[s_name]["hist_read_length"]['X']
            Y = self.sequana_data[s_name]["hist_read_length"]['Y']
            #Y = [y / sum(Y) for y in Y]
            data[s_name] = {x:y for x,y in zip(X[1:], Y)}
            try:
                X = self.sequana_data[s_name]["hist_read_length"]['X']
                Y = self.sequana_data[s_name]["hist_read_length"]['Y']
                #Y = [y / sum(Y) for y in Y]
                data[s_name] = {x:y for x,y in zip(X[1:], Y)}
            except KeyError:
                pass

        if len(data) == 0:
            logger.debug('no data for the read length plots')
            return None

        pconfig = {
            'id': 'sequana_pacbio_hist_length',
            'title': 'Read length',
            'ylab': '#',
            'xlab': 'Length',
            'ymin': 0,
            'xmax': 50000,
            'xmin': 0,
            'yDecimals': False,
            'tt_label': '<b>{point.x}Length</b>: {point.y}',
            #'colors': self.get_status_cols('per_sequence_gc_content'),
            #'data_labels': [
            #    {'name': 'Percentages', 'ylab': 'Percentage'},
            #    {'name': 'length', 'ylab': '#'}
            #]
        }

        self.add_section (
            name = 'Read length histograms',
            anchor = 'fastqc_per_sequence_gc_content',
            description = "Read length histogram",
            plot = linegraph.plot(data, pconfig))

    def parse_logs(self, log_dict):
        import json
        log_dict = json.loads(log_dict)
        data = {}
        data['count'] = log_dict['read_stats']["count"]
        data["mean_length"] = log_dict['read_stats']["mean"]
        data["mean_gc"] = log_dict['mean_gc']
        data["hist_gc"] = log_dict['hist_gc']
        data["hist_read_length"] = log_dict['hist_read_length']
        return data

    def populate_columns(self):
        headers = {}
        if any(['count' in self.sequana_data[s] for s in self.sequana_data]):
            headers['count'] = {
                'title': 'Read count',
                'description': 'read count',
                'min': 0,
                'scale': 'RdYlGn',
                'format': '{:,0d}',
                'shared_key': 'count',
            }

        if any(['mean_gc' in self.sequana_data[s] for s in self.sequana_data]):
            headers["mean_gc"] = {
                'title': 'Mean GC (%)',
                'description': 'mean GC content in %',
                'max': 100,
                'min': 0,
                'scale': 'RdYlGn',
                'format': '{:,.2f}'
            }

        if any(['mean_length' in self.sequana_data[s] for s in self.sequana_data]):
            headers["mean_length"] = {
                'title': 'Mean length',
                'description': 'mean length',
                'max': 15000,
                'min': 0,
                'scale': 'RdYlGn',
                'format': '{:,.2f}'
            }
        if len(headers.keys()):
            self.general_stats_addcols(self.sequana_data, headers)


