# coding: utf-8
#
#  This file is part of Sequana software
#
#  Copyright (c) 2020 - Sequana Development Team
#
#  Distributed under the terms of the 3-clause BSD license.
#  The full license is in the LICENSE file, distributed with this software.
#
#  website: https://github.com/sequana/sequana
#  documentation: http://sequana.readthedocs.io
#
##############################################################################
"""Module to write variant calling report"""
import ast

import pandas as pd

from sequana.modules_report.base_module import SequanaBaseModule
from sequana.utils.datatables_js import DataTable


class Enrichment(SequanaBaseModule):
    """ Write HTML report of variant calling. This class takes a csv file
    generated by sequana_variant_filter.
    """
    def __init__(self, folder):
        """.. rubric:: constructor

        """
        super().__init__()
        self.title = "Enrichment"

        self.create_report_content()
        self.create_html("rnadiff.html")

    def create_report_content(self):
        self.sections = list()

        self.summary()
        self.add_go()
        self.add_kegg()

    def summary(self):
        """ Add information of filter.
        """
        S = self.rnadiff.summary()
        
        self.sections.append({
            'name': "Summary",
            'anchor': 'filters_option',
            'content':
                """
<p>The final Differententially Gene Expression (DGE) analysis
led to {} up and {} down genes (total {})</p>""".format(S.loc['up'][0],
S.loc['down'][0],
S.loc['all'][0])
        })

    def add_kegg(self):
        pass

    def add_go(self):
        style = "width:65%"
        def dendogram(filename):
            import pylab
            pylab.ioff()
            pylab.clf()
            self.rnadiff.plot_dendogram()
            pylab.savefig(filename)
            pylab.close()
        html = """<p>The following image shows a hierarchical
clustering of the whole sample set. The data was log-transformed first.
</p>{}<hr>""".format(
        self.create_embedded_png(dendogram, "filename", style=style))

        self.sections.append({
           "name": "Clusterisation",
           "anchor": "table",
           "content": html 
         })
