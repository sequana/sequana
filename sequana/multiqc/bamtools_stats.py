#!/usr/bin/env python

""" MultiQC module to parse output from sequana"""
import os
import re
from collections import OrderedDict
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
                #name = "bamtools.stats." + name.replace("sequana_bamtools_stats_", "")
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
            log.debug("No samples found: sequana_bamtools_stats")
            raise UserWarning

        log.info("Found {} reports".format(len(self.sequana_data)))

        self.populate_columns()
        self.add_total_read_section()
        self.add_mapped_reads_section()
        self.add_duplicates_section()
        self.add_forward_reverse_section()
        self.add_failed_QC()
        self.add_mapped_vs_unmapped_chart()

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
        # FIXME will fail if instead of 0.001% we have 1e-3%
        data = {}
        regexes = {
            'total_reads': r"Total reads:\s*(\d+)",
            'mapped_reads': r"Mapped reads:\s*(\d+)",
            'mapped_reads_pct': r"Mapped reads:\s*\d+\s+\((.+)%\)",
            'forward_strand': r"Forward strand:\s*(\d+)",
            'forward_strand_pct': r"Forward strand:\s*\d+\s+\((.+)%\)",
            'reverse_strand': r"Reverse strand:\s*(\d+)",
            'reverse_strand_pct': r"Reverse strand:\s*\d+\s+\((.+)%\)",
            'failed_qc': r"Failed QC:\s*(\d+)",
            'failed_qc_pct': r"Failed QC:\s*\d+\s+\((.+)%\)",
            'duplicates': r"Duplicates:\s*(\d+)",
            'duplicates_pct': r"Duplicates:\s*\d+\s+\((.+)%\)",
            'paired_end': r"Paired-end reads:\s*(\d+)",
            'paired_end_pct': r"Paired-end reads:\s*\d+\s+\((.+)%\)",
            'proper_pairs': r"'Proper-pairs'\s*(\d+)",
            'proper_pairs_pct': r"'Proper-pairs'\s*\d+\s+\((.+)%\)",
            'both_mapped': r"Both pairs mapped:\s*(\d+)",
            'both_mapped_pct': r"Both pairs mapped:\s*\d+\s+\((.+)%\)",
            'read_1': r"Read 1:\s*(\d+)",
            'read_2': r"Read 2:\s*(\d+)",
            'singletons': r"Singletons:\s*(\d+)",
            'singletons_pct': r"Singletons:\s*\d+\s+\((.+)%\)",
        }

        import re
        for k, v in regexes.items():
            res = re.findall(v, log_dict)
            if len(res) == 0:
                #logging.warning("Found no entries for {}".format(k))
                data[k] = None
            elif len(res) == 1:
                data[k] = float(res[0])
            else:
                logging.warning("Found several entries for {}".format(k))
                data[k] = None
        return data

    def add_failed_QC(self):
        data = {}
        for name in self.sequana_data.keys():
            data[name] = {'failed_qc_pct': self.sequana_data[name]["failed_qc_pct"]}

        pconfig = {
            "title": "Failed QC (%)",
            "percentages": True,
            "min": 0,
            "max":100,
            "logswitch": False,
        }

        total = sum([data[name]['failed_qc_pct'] for name in self.sequana_data.keys()])
        if total == 0:
            plot = None
            description = "Failed QC (none found)." 
        else:
            plot = bargraph.plot(data, None, pconfig)
            description = "Failed QC (%)." 

        self.add_section(
            name = 'Failed QC (%)',
            anchor = 'failed_qc',
            description = 'Failed QC',
            helptext = "",
            plot = plot)

    def add_forward_reverse_section(self):
        data = {}
        for name in self.sequana_data.keys():
            forward = self.sequana_data[name]["forward_strand_pct"]
            reverse = self.sequana_data[name]["reverse_strand_pct"]
            if reverse is None:
                reverse = 0
            if forward is None:
                forward = 0
            data[name] = {'fwd': forward, "rev":reverse}

        pconfig = {
            "title": "Forward/reverse",
            "percentages": True,
            "min": 0,
            "max":100,
            "logswitch": False,
        }

        keys = OrderedDict()
        keys['fwd'] = {'color': '#437bb1', 'name': 'Forward'}
        keys['rev'] = {'color': '#b1084c', 'name': 'Reverse'}

        self.add_section(
            name = 'Forward/Reverse',
            anchor = 'fwd_rev',
            description = 'Forward reverse',
            helptext = "",
            plot = bargraph.plot(data, keys, pconfig))

    def add_total_read_section(self):
        data = {}
        for name in self.sequana_data.keys():
            data[name] = {'total_read': self.sequana_data[name]["total_reads"]}

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
            mapped_reads = self.sequana_data[name]["mapped_reads_pct"]
            data[name] = {'mapped_reads_pct': mapped_reads}
        pconfig = {
            "title": "Mapped reads (%)",
            "percentages": False,
            "min": 0,
            "max": 100,
            "logswitch": True,
        }

        self.add_section(
            name = 'Mapped reads (%)',
            anchor = 'mapped reads',
            description = 'Total mapped reads found in the BAM file.',
            helptext = "",
            plot = bargraph.plot(data, None, pconfig))

    def add_duplicates_section(self):
        data = {}
        for name in self.sequana_data.keys():
            duplicates_pct = float(self.sequana_data[name]["duplicates_pct"])
            data[name] = {'duplicates_pct': duplicates_pct}

        pconfig = {
            "title": "% Duplicates (%)",
            "description": "% Duplicated Reads",
            "suffix": "%",
            "min": 0,
            "max": 100,
            "scale": "OrRd",
            "logswitch": False,
        }

        total = sum([data[name]['duplicates_pct'] for name in self.sequana_data.keys()])
        if total==0:
            plot = None
            description = "Duplicated reads (none found)." 
        else:
            plot = bargraph.plot(data, None, pconfig)
            description = "Duplicated reads." 

        self.add_section(
            name = 'Duplicates (%)',
            anchor = 'duplicates',
            description = description,
            helptext = "",
            plot = plot)

    def add_mapped_vs_unmapped_chart(self):
        data = {}
        for name in self.sequana_data.keys():
            total = float(self.sequana_data[name]["total_reads"])
            mapped = float(self.sequana_data[name]["mapped_reads"])
            unmapped = total - mapped
            data[name] = {
                'reads_mapped': 100 * mapped/total,
                'reads_unmapped': 100 * unmapped/total
            }

        pconfig = {
            "title": "Mapped vs unmapped reads (%)",
            "percentages": True,
            "min": 0,
            "max": 100,
            "logswitch": False,
        }

        keys = OrderedDict()
        keys['reads_mapped'] = {'color': '#437bb1', 'name': 'Mapped'}
        keys['reads_unmapped'] = {'color': '#b1084c', 'name': 'Unmapped'}

        self.add_section(
            name = "Mapped vs unmapped reads (%)",
            anchor = "bamtools_stats_alignment",
            description  = "Alignment metrics",
            plot = bargraph.plot(data, keys, pconfig)
        )

    def populate_columns(self):
        headers = {}

        if any(['total_reads' in self.sequana_data[s] for s in self.sequana_data]):
            headers['total_reads'] = {
                'title': 'Read count',
                'description': 'read count',
                'min': 0,
                'scale': 'RdYlGn',
                'format': '{0:.0f}',
                'shared_key': 'count',
            }

        if any(['total_reads' in self.sequana_data[s] for s in self.sequana_data]):
            headers['mapped_reads_pct'] = {
                'title': 'Mapped reads (%)',
                'description': 'mapped reads (%)',
                'min': 0,
                'max': 100,
                'scale': 'RdYlGn',
                'format': '{0:.3f}',
                'shared_key': 'count',
            }



        if len(headers.keys()):
            self.general_stats_addcols(self.sequana_data, headers)

