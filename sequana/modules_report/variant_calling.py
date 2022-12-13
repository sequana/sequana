# coding: utf-8
#
#  This file is part of Sequana software
#
#  Copyright (c) 2016 - Sequana Development Team
#
#  File author(s):
#      Dimitri Desvillechabrol <dimitri.desvillechabrol@pasteur.fr>,
#          <d.desvillechabrol@gmail.com>
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

from sequana.lazy import pandas as pd

from sequana.modules_report.base_module import SequanaBaseModule
from sequana.utils.datatables_js import DataTable


class VariantCallingModule(SequanaBaseModule):
    """Write HTML report of variant calling. This class takes a csv file
    generated by sequana_variant_filter.
    """

    def __init__(self, data):
        """.. rubric:: constructor

        :param data: it can be a csv filename created by
        sequana.freebayes_vcf_filter or a
        :class:`freebayes_vcf_filter.Filtered_freebayes` object.
        """
        super().__init__()
        self.title = "Variant Calling Report"
        try:
            with open(data, "r") as fp:
                self.filename = data
                line = fp.readline()
                if line.startswith("# sequana_variant_calling"):
                    string_dict = line.split(";")[-1].strip()
                    try:
                        self.filter_dict = ast.literal_eval(string_dict)
                    except SyntaxError:
                        self.filter_dict = None
                    self.df = pd.read_csv(fp)
        except FileNotFoundError:
            msg = (
                "The csv file is not present. Please, check if your" " file is present."
            )
            raise FileNotFoundError(msg)
        except TypeError:
            self.df = data.df
            self.filter_dict = data.vcf.filters_params
        self.create_report_content()
        self.create_html("variant_calling.html")

    def create_report_content(self):
        self.sections = list()

        if self.filter_dict:
            self.filters_information()
        self.variant_calling()

    def filters_information(self):
        """Add information of filter."""
        self.sections.append(
            {
                "name": "Filter Options",
                "anchor": "filters_option",
                "content": "<p>All filters parameters used is presented in this list:</p>"
                "\n<ul><li>freebayes_score: {freebayes_score}</li>\n"
                "<li>frequency: {frequency}</li>\n"
                "<li>min_depth: {min_depth}</li>\n"
                "<li>forward_depth: {forward_depth}</li>\n"
                "<li>reverse_depth: {reverse_depth}</li>\n"
                "<li>strand_ratio: {strand_ratio}</li></ul>\n"
                "Note:<ul><li>frequency: alternate allele / depth</li>\n"
                "<li>min_depth: minimum alternate allele present</li>\n"
                "<li>forward_depth: minimum alternate allele present on "
                "forward strand</li>\n"
                "<li>reverse_depth: minimum alternate allele present on "
                "reverse strand</li>\n"
                "<li>strand_ratio: alternate allele forward / (alternate "
                "allele forward + alternate allele reverse)</li>"
                "</ul>".format(**self.filter_dict),
            }
        )

    def variant_calling(self):
        """Variants detected section."""
        datatable = DataTable(self.df, "vc")
        # set options
        datatable.datatable.datatable_options = {
            "scrollX": "true",
            "pageLength": 30,
            "scrollCollapse": "true",
            "dom": "Bfrtip",
            "buttons": ["copy", "csv"],
        }
        js = datatable.create_javascript_function()
        html_tab = datatable.create_datatable(float_format="%.3f")
        self.sections.append(
            {
                "name": "Variants Detected",
                "anchor": "basic_stats",
                "content": "<p>This table gives variants detected by freebayes after "
                "filtering. The important metrics are the depth (if not enough reads supports "
                "the variant it should be ignored); the strand_balance (forward and reverse number "
                " of reads supporting the variants should be similar e.g balance of 0.5); the fisher "
                "pvalue (variants with pvalue<0.05 should be rejected since the strand balance of"
                "alternate and reference are different).</p><p>Note: the freebayes score can be"
                " understood as 1 - P(locus is homozygous given the data)</p> {0}\n{1}\n".format(
                    js, html_tab
                ),
            }
        )
