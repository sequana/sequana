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
"""Module to write enrichment report"""
import ast
import os
import sys

from sequana.lazy import pandas as pd
from sequana.lazy import pylab

from sequana.modules_report.base_module import SequanaBaseModule
from sequana.utils.datatables_js import DataTable

from sequana import logger
logger.name = __name__


class Enrichment(SequanaBaseModule):
    """ Write HTML report of variant calling. This class takes a csv file
    generated by sequana_variant_filter.
    """
    def __init__(self, rnadiff_folder, taxon,
                 kegg_organism=None,
                 enrichment_params={
                        "padj": 0.05,
                        "log2_fc": 3,
                        "kegg_background": None,
                        "mapper": None,
                        "preload_directoryr": None,
                        },
                go_only=False,
                kegg_only=False,
                command=""
                ):
        """.. rubric:: constructor

        """
        super().__init__()
        self.title = "Enrichment"

        self.command = command
        self.rnadiff_folder = rnadiff_folder
        self.enrichment_params = enrichment_params
        self.taxon = taxon
        if taxon == 10090:
            self.organism = "mmu"
        elif taxon == 9606:
            self.organism = "hsa"
        else:
            if kegg_organism is None:
                logger.error("You must specify the kegg organism name if not human or mouse: eg., eco for ecoli")
                # figure out the organism from taxon 
                raise NotImplementedError
            else:
                from bioservices import KEGG
                k = KEGG()
                k.organism = kegg_organism # validates the organism name
                self.organism = kegg_organism

        if self.enrichment_params['preload_directory']:
            pathname = self.enrichment_params['preload_directory']
            if os.path.exists(pathname) is False:
                logger.error("{} does not exist".format(pathname))
                sys.exit(1)


        from sequana.rnadiff import RNADiffResults
        self.rnadiff = RNADiffResults(self.rnadiff_folder)

        self.create_report_content(go_only=go_only, kegg_only=kegg_only)
        self.create_html("enrichment.html")

    def create_report_content(self, go_only=False, kegg_only=True):
        self.sections = list()
        self.summary()
        if kegg_only is True:
            self.add_kegg()
        elif go_only is True:
            self.add_go()
        else:
            self.add_go()
            self.add_kegg()
        self.sections.append({
            'name': "4 - Info",
            'anchor': 'command',
            'content': self.command})

    def summary(self):
        """ Add information of filter."""

        S = self.rnadiff.summary()
        Sup = S.loc['up'][0]
        Sdown = S.loc['down'][0]
        Stotal = Sup + Sdown
        try: # if it exists already, do not copy 
            link_rnadiff = self.copy_file(self.rnadiff.filename, ".")
        except:
            link_rnadiff = self.rnadiff.filename
        log2fc = self.enrichment_params["log2_fc"]

        self.sections.append({
            'name': "1 - Summary",
            'anchor': 'filters_option',
            'content':
                f"""
<p>The final Differententially Gene Expression (DGE) analysis
led to {Sup} up and {Sdown} down genes (total {Stotal})</p>

<p>In the following sections, you will find the KEGG Pathway enrichment and GO
terms enrichment. The input data for those analyis is the output of the RNADiff
analysis where adjusted p-values above 0.05 are excluded. Moreover, we removed 
candidates with log2 fold change below {log2fc}. In the following plots you can
find the first GO terms that are enriched, keeping a maximum of 50 identifiers. 
</p>


<p>The input file from the RNADiff analysis is downloadable <a
href="{link_rnadiff}">here</a>.</p>
"""
        })

    def add_go(self):   
        # somehow, logger used here and in add_kegg must be global. If you call
        # add_go and then add_kegg, logger becomes an unbound local variable. 
        # https://stackoverflow.com/questions/10851906/python-3-unboundlocalerror-local-variable-referenced-before-assignment
        global logger
        logger.info("Enrichment module: go term")
        style="width:85%"
        level = logger.level
        logger.level = "INFO"
        from sequana.enrichment import PantherEnrichment
        self.pe = PantherEnrichment(self.rnadiff_folder, self.taxon,
            log2_fc_threshold=self.enrichment_params['log2_fc'])

        # create html table for taxon information
        _taxon_id = self.pe.taxon_info['taxon_id']
        _taxon_name = self.pe.taxon_info['long_name']


        html_intro = f"""<div><p>The taxon used is {_taxon_name} (ID {_taxon_id}).<br>"""
        html_intro += """
        <br>The user input filter excluded adjusted p-values below {} and fold
change in the range [{}, {}] keeping a total of {} down genes and {} genes.</p>
<p>The first plot gather the main molecular function, biological process and
cellular components altogether while the 3 next plots focus on each of those 3
categories. </p>
        </div>
        """.format(
                self.pe.summary['padj_threshold'],
                self.pe.summary['fold_change_range'][0],
                self.pe.summary['fold_change_range'][1],
                self.pe.summary['DGE_after_filtering']['down'],
                self.pe.summary['DGE_after_filtering']['up']
                )

        # compute enrichment. This may take time.
        self.pe.compute_enrichment_down(ontologies=[self.pe.MF, self.pe.BP, self.pe.CC])
        self.pe.compute_enrichment_up(ontologies=[self.pe.MF, self.pe.BP, self.pe.CC])

        # a utility function to create the proper html table
        def get_html_table(this_df, identifier):
            df = this_df.copy()
            links = ["https://www.ebi.ac.uk/QuickGO/term/{}".format(x) for x in df["id"]]
            df['links'] = links
            for x in ['term', 'fdr2', 'abs_log2_fold_enrichment', 'pct_diff_expr']:
                try:del df[x]
                except:pass

            first_col = df.pop("id")
            df.insert(0, "id", first_col)
            df = df.sort_values(by="fold_enrichment", ascending=False)

            datatable = DataTable(pd.DataFrame(df), identifier)
            datatable.datatable.set_links_to_column("links", "id")
            datatable.datatable.datatable_options = {
                 'scrollX': 'true',
                 'pageLength': 10,
                 'scrollCollapse': 'true',
                 'dom': 'Bfrtip',
                 'buttons': ['copy', 'csv']
            }
            js = datatable.create_javascript_function()
            html_table = datatable.create_datatable(float_format='%E')
            return js + html_table

        all_cases = ['GO:0003674', 'GO:0008150', 'GO:0005575']
        case_names = ['MF_BP_CC', 'MF', 'BP', 'CC']

        # all down cases
        self._temp_df = {}
        self._minus = {}
        self._plus = {}

        html_down = ""
        for case, case_name in zip([all_cases, self.pe.MF, self.pe.BP, self.pe.CC], case_names):
            def plot_go_terms_down(filename, ontologies, case_name):
                df = self.pe.plot_go_terms("down", ontologies=ontologies,
                                                  compute_levels=False)
                self._temp_df[case_name] = df.copy()
                self._plus[case_name] = sum(df.plus_minus == '+')
                self._minus[case_name] = sum(df.plus_minus == '-')
                pylab.savefig("Panther_down_{}.png".format(case_name))
                pylab.savefig(filename)
                pylab.close()
            image = self.create_embedded_png(plot_go_terms_down, "filename",
                                             style=style, ontologies=case,
                                             case_name=case_name)
            html_down += "<h4>Down - {}</h4><p>For {} ({}), we found {} go terms. Showing 50 here below (at most). The full list is downlodable from the CSV file hereafter.</p>".format(case_name, case_name, case, 
                self._plus[case_name]+ self._minus[case_name]) 
            html_down += image + " <br>"
            html_down += get_html_table(self._temp_df[case_name], "GO_table_{}".format(case_name))
        
        filenames = []
        for case, case_name in zip([self.pe.MF, self.pe.BP, self.pe.CC], ["MF","BP","CC"]):
            filename = "Chart_down_{}.png".format(case_name)
            if len(self._temp_df[case_name]):
                logger.info("Saving chart for case {} (down)".format(case_name))
                self.pe.save_chart(self._temp_df[case_name], filename)
                filenames.append(filename)
        fotorama_down = "<h4>Charts down</h2>" + self.add_fotorama(filenames)
        
        # all up cases
        self._temp_df = {}
        self._minus = {}
        self._plus = {}
        html_up = ""
        for case, case_name in zip([all_cases, self.pe.MF, self.pe.BP, self.pe.CC], case_names):
            def plot_go_terms_up(filename, ontologies, case_name):
                df = self.pe.plot_go_terms("up", ontologies=ontologies,
                                                 compute_levels=False)
                self._temp_df[case_name] = df.copy()
                self._plus[case_name] = sum(df.plus_minus == '+')
                self._minus[case_name] = sum(df.plus_minus == '-')
                pylab.savefig("Panther_up_{}.png".format(case_name))
                pylab.savefig(filename)
                pylab.close()
            image = self.create_embedded_png(plot_go_terms_up, "filename",
                                             style=style, ontologies=case,
                                             case_name=case_name)
            html_up += "<h4>Up - {}</h4><p>For {} ({}), we found {} go terms. Showing 50 here below (at most). The full list is downlodable from the CSV file hereafter.</p>".format(case_name, case_name, case, 
                self._plus[case_name]+ self._minus[case_name]) 
            html_up += image + " <br>"

        filenames = []
        for case, case_name in zip([self.pe.MF, self.pe.BP, self.pe.CC], ["MF","BP","CC"]):
            filename = "Chart_up_{}.png".format(case_name)
            if len(self._temp_df[case_name]):
                logger.info("Saving chart for case {} (up)".format(case_name))
                self.pe.save_chart(self._temp_df[case_name], filename)
                filenames.append(filename)
        fotorama_up = "<h4>Charts up</h2>" + self.add_fotorama(filenames)

        html = f"{html_intro} {html_down} <hr> {fotorama_down}<hr> {html_up} <hr> {fotorama_up}<hr>"
        self.sections.append({"name": "2 - GO", "anchor": "go", "content": html})
        logger.level = level

    def add_kegg(self):
        global logger
        logger.info("Enrichment module: kegg term")
        style="width:45%"
        from sequana.enrichment import KeggPathwayEnrichment

        ke = KeggPathwayEnrichment(self.rnadiff,
            self.organism,
            mapper=self.enrichment_params["mapper"],
            log2_fc=self.enrichment_params['log2_fc'],
            background=self.enrichment_params['kegg_background'],
            preload_directory=self.enrichment_params['preload_directory'])

        logger.info("Saving all pathways in kegg_pathways/mmu")
        ke.export_pathways_to_json()

        # Image kegg pathways down
        def plot_barplot_down(filename):
            ke.barplot('down')
            pylab.savefig(filename)
        img_barplot_down = self.create_embedded_png(plot_barplot_down, "filename", style=style)
        def plot_scatter_down(filename):
            ke.scatterplot('down')
            pylab.savefig(filename)
        img_scatter_down = self.create_embedded_png(plot_scatter_down, "filename", style=style)

        # Image kegg pathways up
        def plot_barplot_up(filename):
            ke.barplot('up')
            pylab.savefig(filename)
        img_barplot_up = self.create_embedded_png(plot_barplot_up, "filename", style=style)
        def plot_scatter_up(filename):
            ke.scatterplot('up')
            pylab.savefig(filename)
        img_scatter_up = self.create_embedded_png(plot_scatter_up, "filename", style=style)

        # Results down (pathway info)
        html_before_table = """<p>Enrichment pathways summary</p>"""
        df_down = ke.barplot('down')
        links = ["https://www.genome.jp/dbget-bin/www_bget?path:{}".format(x) for x in df_down["pathway_id"]]
        df_down['links'] = links
        df_down = df_down[["pathway_id", "name", "size", "Overlap", "P-value", 
            "Adjusted P-value", "Genes", "links"]]

        # save pathways and add fotorama
        from sequana import logger
        level = logger.level
        logger.level = "WARNING"
        from easydev import Progress
        pb = Progress(len(df_down))
        files = []
        for i, ID in enumerate(df_down['pathway_id']):
            df = ke.save_pathway(ID)
            files.append(ID + ".png")
            pb.animate(i+1)
        fotorama_down = self.add_fotorama(files, width=800)


        datatable = DataTable(df_down, 'kegg_down')
        datatable.datatable.set_links_to_column("links", "pathway_id")
        datatable.datatable.datatable_options = {
             'scrollX': 'true',
             'pageLength': 20,
             'scrollCollapse': 'true',
             'dom': 'Bfrtip',
             'buttons': ['copy', 'csv']
        }
        js_table_down = datatable.create_javascript_function()
        html_table_down = datatable.create_datatable(float_format='%E')


        # Results up (pathway info)
        df_up = ke.barplot('up')
        links = ["https://www.genome.jp/dbget-bin/www_bget?path:{}".format(x) for x in df_up["pathway_id"]]
        df_up['links'] = links
        df_up = df_up[["pathway_id", "name", "size", "Overlap", "P-value", "Adjusted P-value", "Genes", "links"]]
        datatable = DataTable(df_up, 'kegg_up')
        datatable.datatable.set_links_to_column("links", "pathway_id")
        datatable.datatable.datatable_options = {
             'scrollX': 'true',
             'pageLength': 20,
             'scrollCollapse': 'true',
             'dom': 'Bfrtip',
             'buttons': ['copy', 'csv']
        }
        js_table_up = datatable.create_javascript_function()
        html_table_up = datatable.create_datatable(float_format='%E')
        pb = Progress(len(df_up))
        files = []
        for i, ID in enumerate(df_up['pathway_id']):
            df = ke.save_pathway(ID)
            files.append(ID + ".png")
            pb.animate(i+1)
        fotorama_up = self.add_fotorama(files, width=800)
        logger.level = level

        Ndown = len(df_down)
        Nup = len(df_up)

        html = f"""
<h3>2.1 - KEGG pathways down regulated</h3>
<p>{Ndown} KEGG pathways are found to be down regulated</p>
<br>
{img_barplot_down}
{img_scatter_down}
<hr>
{js_table_down} {html_table_down}
<hr>
{fotorama_down}


<h3>2.1 - KEGG pathways up regulated</h3>
<p>{Nup} KEGG pathways are found to be up regulated</p>
<br>
{img_barplot_up}
{img_scatter_up}
<hr>
{js_table_up} {html_table_up}
<hr>
{fotorama_up}
"""
        self.sections.append({"name": "2 - KEGG", "anchor": "kegg", "content": html})



