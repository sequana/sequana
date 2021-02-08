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
log = colorlog.getLogger(__name__)



class MultiqcModule(BaseMultiqcModule):

    def __init__(self):

        # Initialise the parent object
        super(MultiqcModule, self).__init__(
            name='Sequana/isoseq_qc',    # name that appears at the top
            anchor='sequana_isoseq_qc',  # ??
            target='sequana_isoseq_qc',  # Name show that link to the following href
            href='http://github.com/sequana/sequana/',
            info="pipelines multi Summary")

        self.sequana_data = {}
        for myfile in self.find_log_files("sequana_isoseq_qc"):
            name = myfile['s_name']
            try:
                parsed_data = self.parse_logs(myfile["f"])
            except Exception as err:
                print("{} could not be parsed".format(name))
                print(err)
                continue
            name = parsed_data["s_name"]
            self.sequana_data[name] = parsed_data

        if len(self.sequana_data) == 0:
            log.debug("No samples found: sequana_isoseq_qc")
            raise UserWarning

        info = "<ul>"
        for this in sorted(self.sequana_data.keys()):
            info += '<li><a href="{}/summary.html">{}</a></li>'.format(this,this,this)
        info += "</ul>"
        href="http://sequana.readthedocs.io/en/master/"
        target = "Sequana"
        mname = '<a href="{}" target="_blank">{}</a> individual report pages:'.format(href, target)
        self.intro = '<p>{} {}</p>'.format( mname, info)


        log.info("Found {} reports".format(len(self.sequana_data)))

        self.populate_columns()
        self.add_productivity()
        self.add_total_bases()

    def add_total_bases(self):
        data = {}
        for name in self.sequana_data.keys():
            data[name] = {
                'total_bases_gb': self.sequana_data[name]["total_bases_gb"],
            }
        pconfig = {
            "title": "Total bases (Gb)",
            "percentages": False,
            "min": 100,
            "logswitch": True,
        }
        self.add_section(
            name = 'total bases (Gb)',
            anchor = 'total_bases_gb',
            description = '',
            helptext = "",
            plot = bargraph.plot(data, None, pconfig))

    def add_productivity(self):
        from collections import OrderedDict
        data = OrderedDict()
        for name in self.sequana_data.keys():
            data[name] = OrderedDict({
                'P0': self.sequana_data[name]["P0"],
                'P1': self.sequana_data[name]["P1"],
                'P2': self.sequana_data[name]["P2"]
            })

        pconfig = {
            "title": "Productivity",
            "percentages": True,
            "min": 0,
            "max": 0,
            "logswitch": False,
        }
        self.add_section(
            name = 'Productivity',
            anchor = 'productivity',
            description = 'P0/P1/P2 productivities',
            helptext = "",
            plot = bargraph.plot(data, None, pconfig))

    def parse_logs(self, log_dict):
        import json
        log_dict = json.loads(log_dict)
        data = {}
        data["s_name"] = log_dict['sample_name']

        log_dict = log_dict["data"]["QC"]
        #data['count'] = log_dict['data']["count"]
        data["P0"] = log_dict["productivity"]["P0"]
        data["P1"] = log_dict["productivity"]["P1"]
        data["P2"] = log_dict["productivity"]["P2"]
        data["total_bases_gb"] = log_dict["total_bases_gb"]
        data["read_length_poly_AVG"] = log_dict["read_length"]["polymerase"]["AVG"]
        data["read_length_poly_N50"] = log_dict["read_length"]["polymerase"]["N50"]

        return data

    def populate_columns(self):
        from collections import OrderedDict
        headers = OrderedDict()

        if any(["total_bases_gb" in self.sequana_data[s] for s in self.sequana_data]):
            headers["total_bases_gb"] = {
                    'title': "total bases (Gb)",
                    'description': "Total bases (Gb)",
                    'min': 0,
                    'max': 100,
                    'scale': 'RdYlGn',
                    'format': '{:,.0d}',
                    'shared_key': 'productivity',
                }

        for this in ["P0", "P1", "P2"]:
            if any([this in self.sequana_data[s] for s in self.sequana_data]):
               headers[this] = {
                    'title': this,
                    'description': this,
                    'min': 0,
                    'max': 100,
                    'scale': 'RdYlGn',
                    'format': '{:,.0d}',
                    'shared_key': 'productivity',
                }


        if len(headers.keys()):
            self.general_stats_addcols(self.sequana_data, headers)


