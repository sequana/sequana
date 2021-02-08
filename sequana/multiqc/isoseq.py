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
            name='Sequana/isoseq',    # name that appears at the top
            anchor='sequana_isoseq',  # ??
            target='sequana_isoseq',  # Name show that link to the following href
            href='http://github.com/sequana/sequana/',
            info="pipelines multi Summary")

        self.sequana_data = {}
        for myfile in self.find_log_files("sequana_isoseq"):
            name = myfile['s_name']

            try:
                parsed_data = self.parse_logs(myfile["f"])
            except Exception as err:
                print("{} {} could not be parsed".format(name, myfile['f']))
                print(err)
                continue
            name = parsed_data["s_name"]
            self.sequana_data[name] = parsed_data

        if len(self.sequana_data) == 0:
            log.debug("No samples found: sequana_isoseq")
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
        self.add_ccs_reads_section()
        self.add_ccs_mean_length_section()
        self.add_isoforms("hq")
        self.add_isoforms("lq")
        self.add_polyA()
        self.add_primes("three")
        self.add_primes("five")
        self.add_flnc()
        #self.add_ratio_flnc_hq()

    def add_flnc(self):
        data = {}
        for name in self.sequana_data.keys():
            data[name] = {
                'flnc': self.sequana_data[name]["flnc"]
            }
        pconfig = {
            "title": "flnc",
            "logswitch": False,
        }
        self.add_section(
            name = 'Number of flnc',
            anchor = "flnc",
            description = "flnc",
            helptext="",
            plot = bargraph.plot(data, None, pconfig))

    def add_primes(self, label):
        name2 = "{} primes".format(label)
        data = {}
        for name in self.sequana_data.keys():
            data[name] = {
                '3 {}'.format(label):
                    self.sequana_data[name]["{}_primes".format(label)]
            }
        pconfig = {
            "title": name2,
            "logswitch": True,
        }
        self.add_section(
            name = 'Number of {}'.format(name2),
            anchor = name2,
            description = name2,
            helptext="",
            plot = bargraph.plot(data, None, pconfig))

    def add_polyA(self):
        data = {}
        for name in self.sequana_data.keys():
            data[name] = {
                'polyA':
                    self.sequana_data[name]["polyA"]
            }
        pconfig = {
            "title": "polyA",
            "logswitch": True,
        }
        self.add_section(
            name = 'Number of polyA',
            anchor = 'poylA',
            description = 'polyA',
            helptext="",
            plot = bargraph.plot(data, None, pconfig))

    def add_isoforms(self, mode):
        data = {}
        for name in self.sequana_data.keys():
            data[name] = {
                'number_{}_isoforms'.format(mode):
                    self.sequana_data[name]["number_{}_isoforms".format(mode)]
            }
        pconfig = {
            "title": "Number of {} isoforms".format(mode.upper()),
            "logswitch": True,
        }
        self.add_section(
            name = 'Number of {} isoforms'.format(mode.upper()),
            anchor = 'number_{}_isoforms'.format(mode),
            description = 'Number of {} isoforms'.format(mode.upper()),
            helptext="",
            plot = bargraph.plot(data, None, pconfig))

    def add_ccs_reads_section(self):
        data = {}
        for name in self.sequana_data.keys():
            data[name] = {
                'number_ccs_reads': self.sequana_data[name]["number_ccs_reads"]
            }

        pconfig = {
            "title": "Number of CCS reads",
            "percentages": False,
            "min": 100,
            "logswitch": True,
        }
        self.add_section(
            name = 'Number of CCS reads',
            anchor = 'number_ccs_reads',
            description = 'Number of CCS reads',
            helptext = "",
            plot = bargraph.plot(data, None, pconfig))

    def add_ccs_mean_length_section(self):
        data = {}
        for name in self.sequana_data.keys():
            data[name] = {
                'mean': self.sequana_data[name]["mean_length"],
            }
        pconfig = {
                "title": "Mean CCS read length",
            "percentages": False,
            "min": 100,
            "logswitch": True,
        }

        self.add_section(
            name = 'Mean CCS read length',
            anchor = 'mean_ccs_read_length',
            description = 'Mean CCS length of the reads',
            helptext = "",
            plot = bargraph.plot(data, None, pconfig))

    def parse_logs(self, log_dict):
        import json
        log_dict = json.loads(log_dict)
        data = {}
        #data['count'] = log_dict['data']["count"]
        data["mean_length"] = log_dict['data']['CCS']["mean_length"]
        data["number_ccs_reads"] = log_dict['data']['CCS']["number_ccs_reads"]
        data["number_hq_isoforms"] = log_dict['data']['hq_isoform']["N"]
        data["number_lq_isoforms"] = log_dict['data']['lq_isoform']["N"]
        data["flnc"] = log_dict['data']['classification']["full_length_non_chimeric_reads"]
        data["polyA"] = log_dict['data']['classification']["polyA_reads"]
        data["three_primes"] = log_dict['data']['classification']["three_prime_reads"]
        data["five_primes"] = log_dict['data']['classification']["five_prime_reads"]
        #data["ratio_flnc_hq"] = data["flnc"] / data["number_hq_isoforms"]

        data["alldata"] = log_dict
        data["s_name"] = log_dict['sample_name']

        return data

    def populate_columns(self):
        headers = {}
        if any(['mean_length' in self.sequana_data[s] for s in self.sequana_data]):
            headers['mean_length'] = {
                'title': 'CCS mean read length',
                'description': 'CCS mean read length',
                'min': 0,
                'scale': 'RdYlGn',
                'format': '{:,.0d}',
                'shared_key': 'count',
            }

        if any(['number_ccs_reads' in self.sequana_data[s] for s in self.sequana_data]):
            headers["number_ccs_reads"] = {
                'title': 'CCS reads',
                'description': 'Number of CCS reads',
                'min': 0,
                'scale': 'RdYlGn',
                'format': '{:,.0d}'
            }

        for this in ['number_hq_isoforms', 'number_lq_isoforms']:
            if any([this in self.sequana_data[s] for s in self.sequana_data]):
                headers[this] = {
                    'title': " ".join(this.split()),
                    'description': " ".join(this.split()),
                    'min': 0,
                    'scale': 'RdYlGn',
                    'format': '{:,.0d}'
                }
        if any(["polyA" in self.sequana_data[s] for s in self.sequana_data]):
                headers["polyA"] = {
                    'title': "polyA",
                    'description': "",
                    'min': 0,
                    'scale': 'RdYlGn',
                    'format': '{:,.0d}'
                }

        if any(["three_primes" in self.sequana_data[s] for s in self.sequana_data]):
                headers["three_primes"] = {
                    'title': "three prime",
                    'description': "",
                    'min': 0,
                    'scale': 'RdYlGn',
                    'format': '{:,.0d}'
            } 
        if any(["five_primes" in self.sequana_data[s] for s in self.sequana_data]):
                headers["five_primes"] = {
                    'title': "five prime",
                    'description': "",
                    'min': 0,
                    'scale': 'RdYlGn',
                    'format': '{:,.0d}'
            }
        if any(["flnc" in self.sequana_data[s] for s in self.sequana_data]):
                headers["flnc"] = {
                    'title': "flnc",
                    'description': "",
                    'min': 0,
                    'scale': 'RdYlGn',
                    'format': '{:,.0d}'
            }
        if len(headers.keys()):
            self.general_stats_addcols(self.sequana_data, headers)


