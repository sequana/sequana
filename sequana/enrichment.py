# -*- coding: utf-8 -*-
#
#  This file is part of Sequana software
#
#  Copyright (c) 2016-2020 - Sequana Development Team
#
#  File author(s):
#      Thomas Cokelaer <thomas.cokelaer@pasteur.fr>
#      Etienne Kornobis <etienne.kornobis@pasteur.fr>
#
#  Distributed under the terms of the 3-clause BSD license.
#  The full license is in the LICENSE file, distributed with this software.
#
#  website: https://github.com/sequana/sequana
#  documentation: http://sequana.readthedocs.io
#
##############################################################################

from pathlib import Path
import re
import os

from sequana.lazy import pandas as pd
from sequana.lazy import pylab
from sequana.lazy import numpy as np
from sequana import logger
from matplotlib_venn import venn2_unweighted, venn3_unweighted
from sequana.rnadiff import RNADiffResults

try:
    import gseapy
except:
    pass


__all__ = ["KeggPathwayEnrichment"]

class PantherEnrichment():  # pragma: no cover
    """
        p = Panther()
        gene_list = pe.rnadiff.df.query("log2FoldChange>1 and padj<0.05")
        gg = ",".join(gene_list.Name)
        res = p.get_enrichment(gg, 83333, 'ANNOT_TYPE_ID_PANTHER_GO_SLIM_MF')
        res['result']
    """
    def __init__(self):
        pass


class KeggPathwayEnrichment():  # pragma: no cover
    """DRAFT IN PROGRESS

    ::

        pe = PathwayEnrichment("rnadiff", "eco")
        pe.compute_enrichment()
        pe.barplot(pe.enrichment['down'])


        # Save all deregulated pathways found by the enrichment:
        up = pe.save_significant_pathways("up")
        down = pe.save_significant_pathways("down")


    """
    def __init__(self, folder, organism, comparison=None, alpha=0.05, fc=0):
        print("DRAFT in progress")
        from bioservices import KEGG
        self.kegg = KEGG(cache=True)
        self.kegg.organism = organism

        self.read_rnadiff(folder, alpha=alpha, fc=fc)
        choices = list(self.rnadiff.dr_gene_lists.keys())
        if comparison:
            assert comparison in choices
            self.comparison = comparison
        elif len(choices) == 1:
            logger.info("One comparison found and set automatically")
            self.comparison = choices[0]

        self.background = len(self.kegg.list(self.kegg.organism).split("\n"))
        logger.info("Set number of genes to {}".format(self.background))

        self._load_pathways()

    def read_rnadiff(self, folder, alpha=0.05, fc=0):
        self.rnadiff = RNADiffResults(folder ,alpha=alpha, fc=fc)

    def _load_pathways(self):
        # This is just loading all pathways once for all
        logger.info("loading all pathways from KEGG. may take time the first time")
        self.pathways = {}
        from easydev import Progress
        pb = Progress(len(self.kegg.pathwayIds))
        for i, ID in enumerate(self.kegg.pathwayIds):
            self.pathways[ID.replace("path:", "")] = self.kegg.parse(self.kegg.get(ID))
            pb.animate(i+1)

        # Some cleanup
        for ID in self.pathways.keys():
            name = self.pathways[ID]['NAME'][0]
            self.pathways[ID]['NAME'] = name.split(" - ", 1)[0]

        # save gene sets
        self.gene_sets = {}
        for ID in self.pathways.keys():
            res =  self.pathways[ID]
            if "GENE" in res.keys():
                results = [res['GENE'][g].split(';')[0] for g in res['GENE'].keys()]
                self.gene_sets[ID] = results
            else:
                print("SKIPPED (no genes) {}: {}".format(ID, res['NAME']))

        # save all pathways info
        self.df_pathways = pd.DataFrame(self.pathways).T
        del self.df_pathways["ENTRY"]
        del self.df_pathways["REFERENCE"]
        go = [x['GO'] if isinstance(x, dict) and 'GO' in x.keys() else None for x in self.df_pathways.DBLINKS]
        self.df_pathways['GO'] = go
        del self.df_pathways["DBLINKS"]

    def plot_genesets_hist(self, bins=20):
        N = len(self.gene_sets.keys())
        pylab.clf()
        pylab.hist([len(v) for k,v in self.gene_sets.items()],bins=bins, lw=1,
            ec="k")
        pylab.title("{} gene sets".format(N))
        pylab.xlabel("Gene set sizes")
        pylab.grid(True)
        a,b = pylab.xlim()
        pylab.xlim([0, b])

    def compute_enrichment(self, background=None):
        if background is None:
            background = self.background
        self.enrichment = {}
        self.enrichment['up'] = self._enrichr("up", background=background)
        self.enrichment['down'] = self._enrichr("down", background=background)
        self.enrichment['all'] = self._enrichr("all", background=background)

    def _enrichr(self, category, background=None, verbose=True):

        if background is None:
            background = self.background

        if isinstance(category, list):
            gene_list = category
        else:
            assert category in ['up', 'down', 'all']
            gene_list = list(self.rnadiff.dr_gene_lists[self.comparison][category])

        enr = gseapy.enrichr(
            gene_list=gene_list,
            gene_sets=self.gene_sets,
            verbose=verbose,
            background=background,
            outdir="test", no_plot=True)

        return enr

    def _get_final_df(self, df, cutoff=0.05, nmax=10):
        df = df.copy()
        df['name'] = [self.pathways[x]['NAME'] for x in df.Term]
        df['size'] = [len(x.split(";")) for x in df.Genes]
        df = df.sort_values("Adjusted P-value")
        df.reset_index(drop=True, inplace=True)
        df = df[df["Adjusted P-value"]<=cutoff]

        if len(df)<nmax:
            nmax = len(df)
        df = df.iloc[0:nmax]
        df = df.sort_values("Adjusted P-value", ascending=False)
        return df

    def barplot(self, enrich, cutoff=0.05, nmax=10):
        df = self._get_final_df(enrich.results, cutoff=cutoff, nmax=nmax)

        pylab.clf()
        pylab.barh(range(len(df)), -pylab.log10(df['Adjusted P-value']))
        pylab.yticks(range(len(df)), df.name)
        pylab.axvline(1.3, lw=2, ls="--", color="r")
        pylab.grid(True)
        pylab.xlabel("Adjusted p-value (log10)")
        pylab.ylabel("Gene sets")
        a,b = pylab.xlim()
        pylab.xlim([0, b])
        pylab.tight_layout()
        return df

    def scatterplot(self, enrich, cutoff=0.05, nmax=10, gene_set_size=[]):
        df = self._get_final_df(enrich.results, cutoff=cutoff, nmax=nmax)

        pylab.clf()
        pylab.scatter(-pylab.log10(df['Adjusted P-value']), range(len(df)), s=10*df['size'],    
            c=df['size'])
 
        pylab.xlabel("Odd ratio")
        pylab.ylabel("Gene sets")
        pylab.yticks(range(len(df)), df.name)
        a,b = pylab.xlim()
        pylab.xlim([0, b])
        pylab.grid(True)
        ax = pylab.gca()

        M = max(df['size'])
        if M >100:
            l1,l2,l3  = "10", "100", str(M)
        else:
            l1,l2,l3 = str(round(M/3)), str(round(M*2/3)), str(M)

        handles = [
            pylab.Line2D([0], [0],marker="o", markersize=5, label=l1, ls=""),
            pylab.Line2D([0], [0],marker="o", markersize=10, label=l2, ls=""),
            pylab.Line2D([0], [0],marker="o", markersize=15, label=l3, ls="")]
        ax.legend(handles=handles, loc="upper left", title="gene-set size")

        pylab.axvline(1.3, lw=2, ls="--", color="r")
        pylab.tight_layout()
        ax = pylab.colorbar(pylab.gci())
        return df

    def _get_summary_pathway(self, pathway_ID):
        genes = self.df_pathways.loc[pathway_ID]['GENE']
        df_down = self.rnadiff.df.query("padj<=0.05 and log2FoldChange<0")
        df_up = self.rnadiff.df.query("padj<=0.05 and log2FoldChange>0")

        #f_down = self.rnadiff.dr_gene_lists[self.comparison]

        logger.info("Total down-regulated: {}".format(len(df_down)))
        logger.info("Total up-regulated: {}".format(len(df_up)))

        mapper = {}
        for k,v in genes.items():
            mapper[v.split(";")[0]] = k
        self.genes = genes
        self.mapper = mapper
        self.df_down = df_down
        self.df_up = df_up
        summary_names = []
        summary_keggids = []
        summary_types = []
        summary_pvalues = []
        summary_fcs = []
        for name, kegg_id in mapper.items():
            summary_names.append(name)
            summary_keggids.append(kegg_id)

            if name.lower() in [x.lower() for x in df_down.Name]:
                pvalue = -pylab.log10(df_down.query("Name==@name").pvalue.values[0])
                fc = df_down.query("Name==@name").log2FoldChange.values[0]
                summary_fcs.append(fc)
                summary_pvalues.append(pvalue)
                summary_types.append("-")
            elif name.lower() in [x.lower() for x in df_up.Name]:
                pvalue = -pylab.log10(df_up.query("Name==@name").pvalue.values[0])
                summary_pvalues.append(pvalue)
                fc = df_up.query("Name==@name").log2FoldChange.values[0]
                summary_fcs.append(fc)
                summary_types.append("+")
            else:
                summary_pvalues.append(None)
                summary_fcs.append(None)
                summary_types.append("=")

        summary = pd.DataFrame({
                    "type":summary_types,
                    "name": summary_names,
                    "pvalue": summary_pvalues,
                    "fc": summary_fcs,
                    "keggid": summary_keggids
                    })
        summary['description'] = [self.pathways[pathway_ID]['GENE'][x] for x in summary.keggid]
        return summary 

    def _get_colors(self, summary):
        colors = {}
        for index, row in summary.iterrows():
            pvalue = row['pvalue']
            type_ = row['type']
            kegg_id = row['keggid']
            if type_ == "-":
                if pvalue>0 and pvalue<5:
                    colors[kegg_id] = "#FF8C00,black"
                elif pvalue<10:
                    colors[kegg_id] = "#FF0000,black"
                else:
                    colors[kegg_id] = "#B22222%2Cblack"
            elif type_ == "+":
                if pvalue>0 and pvalue<5:
                    colors[kegg_id] = "#9ACD32,black"
                elif pvalue<10:
                    colors[kegg_id] = "#008000,black"
                else:
                    colors[kegg_id] = "#006400,#000000"
            else:
                colors[kegg_id] = "grey,black"
        return colors

    def save_pathway(self, pathway_ID, scale=None, show=False, filename=None):

        summary = self._get_summary_pathway(pathway_ID)
        colors = self._get_colors(summary)

        logger.info("pathway {} total genes: {}".format(pathway_ID, len(summary)))
        count_up = len(summary.query("type == '+'"))
        count_down = len(summary.query("type == '-'"))
        logger.info("this pathway down-regulared genes: {}".format(count_down))
        logger.info("this pathway up-regulated genes: {}".format(count_up))

        url = "https://www.kegg.jp/kegg-bin/show_pathway"
        #dcolor = "white"  --> does not work with the post requests unlike get
        # requests
        params = {"map": pathway_ID,
                "multi_query":  "\r\n".join(["{} {}".format(k,v) for k,v in colors.items()])}

        self.params = params
        import requests
        html_page = requests.post(url, data=params)

        self.tmp = html_page
        html_page = html_page.content.decode()

        links_to_png = [x for x in html_page.split() if "png" in x and x.startswith("src")]
        link_to_png =  links_to_png[0].replace("src=", "").replace('"', '')
        r = requests.get("https://www.kegg.jp/{}".format(link_to_png))

        if filename is None:
            filename = "{}.png".format(pathway_ID)

        with open(filename, "wb") as fout:
            fout.write(r.content)

        return summary

    def save_all_pathways(self):
        # This does not do any enrichment. Just save all pathways once for all
        # with useful information
        for ID in self.kegg.pathwayIds:
            self.save_pathway(ID)

    def save_significant_pathways(self, mode, cutoff=0.05, nmax=20,
        background=None):

        if background is None:
            background = self.background

        # select the relevant pathways
        df = self._enrichr(mode, self.background).results
        df = self._get_final_df(df, cutoff=cutoff, nmax=nmax)
        if len(df) == nmax:
            logger.warning("Restricted pathways to {}".format(nmax))

        logger.info("saving {} deregulated pathways".format(len(df)))

        summaries = {}
        # save them
        for ID in df['Term']:
            summary = self.save_pathway(ID, filename="{}_{}.png".format(ID, mode))
            summaries[ID] = summary
        return summaries


