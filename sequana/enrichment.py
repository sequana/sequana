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


__all__ = ["PantherEnrichment", "KeggPathwayEnrichment"]


class PantherEnrichment():  
    """

        # This will read your rnadiff results and tstore the rlevants genes into 
        # mygenes_u, mygenes_down, mygenes attributes.

        By default, we keep only the genes with a adjusted pvalue <= 0.05. The
        fold_change threshold is on a log2 scale and set to 0 (meaning no
        filtering). Only one value is requested and
        used to filter out positive and negative fold change (on the log2
        scale). In other word, a log2 fold change threshold of 2 means that we 
        filter out values between -2 and 2.

        If you prefer to work in natural scale, just set the parameter
        fc_threshold, which overwrite the log2_fc_threshold parameter.

        ::

            e = PantherEnrichment("rnadiff_folder"), log2_fc_threshold=1)
            e.compute_enrichment(pe.mygenes_down, 83333)

        e.stats contains some statistics. One important is the list of unmapped
        genes. The results from the GO enrichment are stored in the attributes
        enrichment. There, we have again adjusted p-value and a fold enrichment,
        which can in turn be filtered or not.

        You can retrieve the cleaned data using the get_data method.

        You can also plot the GO terms that are significantly enriched using::

            e.plot_go_terms(['GO:0003674', 'GO:0008150', 'GO:0005575'])

        This function returns the dataframe used during the plotting.

        If you want to look at the up regulated genes only,

            e.compute_enrichment(pe.mygenes_up, 83333)
            e.plot_go_terms(['GO:0003674', 'GO:0008150', 'GO:0005575'])


        df = e.plot_go_terms(['GO:0003674', 'GO:0008150', 'GO:0005575'],
                log=False, include_negative_enrichment=False,
                fontsize=8, sort_by='fold_enrichment',
                show_pvalues=True, fdr_threshold=0.05)


    """
    def __init__(self, folder, requests_per_sec=10, padj_threshold=0.05,
        log2_fc_threshold=0, fc_threshold=None,
        pattern="*complete*.xls"):

        assert log2_fc_threshold>=0, "log2 fc_threshold must be >=0"

        if fc_threshold is not None:
            log2_fc_threshold = pylab.log2(fc_threshold)

        from bioservices import panther, quickgo, uniprot
        self.panther = panther.Panther()
        self.valid_taxons = [x['taxon_id'] for x in self.panther.get_supported_genomes()]

        self.quickgo = quickgo.QuickGO(cache=True)
        self.uniprot = uniprot.UniProt(cache=True)

        self.quickgo.requests_per_sec = requests_per_sec
        self.uniprot.requests_per_sec = requests_per_sec

        self.ancestors = {"MF": "GO:0003674", "CC": "GO:0005575", "BP": "GO:0008150"}
        self.aspects = {"MF": "molecular_function"}
        self.ontologies = [
            'GO:0003674', 'GO:0008150', 'GO:0005575',
            'ANNOT_TYPE_ID_PANTHER_GO_SLIM_MF',
            'ANNOT_TYPE_ID_PANTHER_GO_SLIM_BP',
            'ANNOT_TYPE_ID_PANTHER_GO_SLIM_CC',
            'ANNOT_TYPE_ID_PANTHER_PC',
            'ANNOT_TYPE_ID_PANTHER_PATHWAY',
            'ANNOT_TYPE_ID_REACTOME_PATHWAY']

        self.ontology_aliases = ["MF", "BP", "CC",
            'SLIM_MF', 'SLIM_BP', 'SLIM_CC', 'PROTEIN', 'PANTHER_PATHWAY',
            'REACTOME_PATHWAY']

        from sequana.rnadiff import RNADiffResults
        self.rnadiff = RNADiffResults(folder, pattern=pattern)
        logger.info("Ignoring pvalue adjusted > {} and fold change in [{}, {}]".format(
            padj_threshold, 1/(2**log2_fc_threshold), 2**log2_fc_threshold ))

        fc_threshold = log2_fc_threshold

        self.mygenes = self.rnadiff.df.query(
            "padj<=@padj_threshold and (log2FoldChange<=-@fc_threshold or log2FoldChange>=@fc_threshold)")
        self.mygenes_down = self.rnadiff.df.query(
            "padj<=@padj_threshold and log2FoldChange<=-@fc_threshold")
        self.mygenes_up = self.rnadiff.df.query(
            "padj<=@padj_threshold and log2FoldChange>=@fc_threshold")

        self.mygenes = list(self.mygenes.sort_values('padj').index)
        self.mygenes_down = list(self.mygenes_down.sort_values('padj').index)
        self.mygenes_up = list(self.mygenes_up.sort_values('padj').index)

        # When using ENSEMBL, prefix "gene:"  should be removed to be understood
        # by PantherDB
        self.mygenes = [x.replace("gene:", "") for x in self.mygenes]
        self.mygenes_down = [x.replace("gene:", "") for x in self.mygenes_down]
        self.mygenes_up = [x.replace("gene:", "") for x in self.mygenes_up]

        logger.info("Kept {} genes ({} up; {} down)".format(
            len(self.mygenes),
            len(self.mygenes_down),
            len(self.mygenes_up)))

    def compute_enrichment(self, mygenes, taxid, ontologies=None, 
            enrichment_test="FISHER", correction="FDR"):
        #taxid=83333 # ecoli
        """
        :param enrichment_test: Fisher or Binomial
        :param correction: FDR or Bonferonni

        The field **number_in_reference** indicates from the reference, the number
        of genes that have a given ontolgy term. For instance, 998 genes have
        the term. This is stored in **number_in_reference**. If the reference 
        contains 4391 genes, and you provided 49
        genes , the **expected** number of genes that have this ontology term is
        49*998/4391 that is 11.1369, which is stored in **"expected**. 

        Now, if you actually find that 14 out of 49
        genes have the term, you need to compare the numbers 11.1369 and 14. Are
        they really different ? The ratio 14 / 11.1369 is stored in
        **fold_enrichment**. The pvalue and FDR are stored as well. 


        Some genes may be missing If you provide 50 genes, you may end up with
        only 45 being mapped onto panther db database. This may explain some
        differences with the expected value.

        Fold enrichment is the number_in_list / expected ratio. Another close metric is the
        fractional difference: (observed - expected) / expected. This metric is
        slighlty less than the fold enrichment

        To get the number f genes from the database, simply take one output, and
        compute number_in_reference * number_in_list / expected


        The fold enrichment is also called odd-ratio.
        """
        if isinstance(mygenes, list):
            mygenes = ",".join(mygenes)

        if ontologies is None:
            ontologies = self.ontologies
        else:
            for x in ontologies:
                assert x in self.ontologies

        self.enrichment = {}
        self.N_genes = len(mygenes)
        for ontology in ontologies:
            logger.info("Compute enrichment for {}".format(ontology))
            results = self.panther.get_enrichment(mygenes, taxid, ontology,
                enrichment_test=enrichment_test, correction=correction)
            #results = [x for x in results['result'] if x['term']['label'] != "UNCLASSIFIED"]

            if isinstance(results["result"], dict):  # pragma: no cover
                results["result"] = [results["result"]]
            pvalues = [x['pValue'] for x in results['result']]
            import statsmodels
            import statsmodels.stats.multitest
            if correction == "FDR":
                fdr = statsmodels.stats.multitest.multipletests(pvalues,
                    method='fdr_bh')[1]
            elif correction.lower() == "bonferroni":
                fdr = statsmodels.stats.multitest.multipletests(pvalues,
                    method='bonferroni')[1]
            for i, pvalue in enumerate(pvalues):
                results["result"][i]['fdr2'] = fdr[i]
                if enrichment_test.lower() == "binomial":
                    results["result"][i]['fdr'] = fdr[i]

            self.enrichment[ontology] = results
        self.stats = dict([(k,len(v['result'])) for k,v in self.enrichment.items()])
        self.stats['input_genes'] = len(mygenes.split(','))
        unmapped = self.enrichment[ontologies[0]]["input_list"]['unmapped_id']
        self.stats['unmapped_genes'] = unmapped
        self.stats['N_unmapped_genes'] = len(unmapped)


        # Here, looking at the FDr, it appears that when using bonferroni,
        # all FDR are set to zeros. Moreover, when using Fisher tests and
        # FDR (supposibly a FDR_BH, the results are noisy as compare to a
        # test from statsmodels. Moreover, when using binomial test, the FDR
        # is not computed... So, we will recompute the FDR ourself

    def get_functional_classification(self, mygenes, taxon): #pragma: no cover ; too slow
        if isinstance(mygenes, list):
            mygenes = ",".join(mygenes)

        res = self.panther.get_mapping(mygenes, taxon)
        res = res["mapped"]
        for i, item in enumerate(res):
            accession = item['accession']
            res[i]['persistent_id'] = self._get_name_given_accession(accession)
        return res

    def _get_name_given_accession(self, accession):  #pragma: no cover
        from bioservices import UniProt
        acc = [x for x in accession.split("|") if x.startswith("UniProtKB")]
        acc = acc[0].split("=")[1]
        res = self.uniprot.get_df(acc, limit=1)
        name = res['Gene names  (primary )'][0]
        return name

    def __get_goinfo(self, mygenes, taxon, annotation_list=["GO:0003674"]): #pragma: no cover

        if isinstance(annotation_list, str):
            annotation_list = [annotation_list]
        assert isinstance(annotation_list, list)

        info = self.get_functional_classification(mygenes, taxon)

        from collections import defaultdict
        goids = defaultdict(list)
        for x in info:
            name = x['persistent_id']
            if "annotation_type_list" in x.keys(): 
                data = x["annotation_type_list"]["annotation_data_type"] 
                for annot in data: 
                    if 'content' not in annot:
                        print(x, annot)
                    print(annot['content'])
                    if annot["content"] in annotation_list: 
                        if "annotation" in annot['annotation_list']: 
                            subdata = annot["annotation_list"]["annotation"] 
 
                            if isinstance(subdata,dict): 
                                goids[name].append(subdata["id"]) 
                            else: 
                                for item in subdata: 
                                    goids[name].append(item["id"]) 
        return goids

    def __get_names(self, mapping_results): #pragma: no cover
        """


            e = PantherEnrichment()
            res = e.panther.get_mapping(mygenes, taxon)
            e.get_names(res['mapped'])

        """
        from bioservices import UniProt
        names = []
        accessions = [x["accession"] for x in mapping_results]
        for accession in accessions:
            acc = [x for x in accession.split("|") if x.startswith("UniProtKB")]
            acc = acc[0].split("=")[1]
            res = self.uniprot.get_df(acc, limit=1)
            name = res['Gene names  (primary )'][0]
            names.append(name)
        return names

    def plot_piechart(self, df):
        # Here we show the GO terms that have number in list > 0
        # Note, that this is dangerous to look only at this picture without
        # the reference plot, which data is not available thourg the pathner API 
        labels = []
        for this in df.query("number_in_list!=0").label.values:
            if len(this)>50:
                labels.append(this)
            else:
                labels.append(this[0:50]+"...")
        pylab.pie(
            df.query("number_in_list!=0").number_in_list,
            labels=labels)
        pylab.tight_layout()

    def get_data(self, ontologies, include_negative_enrichment=True, fdr=0.05):

        if isinstance(ontologies, str):
            ontologies = [ontologies]
        else:
            assert isinstance(ontologies, list)
        # First, we select the required ontologies and build a common data set
        all_data = []
        for ontology in ontologies:
            data = self.enrichment[ontology]['result']
            if isinstance(data, dict):
                # there was only one hit, we expect:
                data = [data]
            all_data.extend(data)
        data = all_data

        # remove unclassified GO terms
        unclassified = [x for x in data if x['term']['label'] == "UNCLASSIFIED"]
        logger.info("Found {} unclassified".format(len(unclassified)))
        data = [x for x in data if x['term']['label'] != "UNCLASSIFIED"]

        df = pd.DataFrame(data)
        if len(df) == 0:
            return df
        else:
            logger.info("Found {} GO terms".format(len(df)))


        df = df.query("number_in_list!=0").copy()
        logger.info("Found {} GO terms with at least 1 gene in reference".format(len(df)))

        # extract the ID and label
        df['id'] = [x['id'] for x in df['term']]
        df['label'] = [x['label'] for x in df['term']]

        # some extra information for convenience
        df["pct_diff_expr"] = df['number_in_list'] *100 / df['number_in_reference']
        df["log2_fold_enrichment"] = pylab.log2(df['fold_enrichment'])
        df["abs_log2_fold_enrichment"] = abs(pylab.log2(df['fold_enrichment']))

        # Some user may want to include GO terms with fold enrichment
        # significanyly below 1 or not. 
        if include_negative_enrichment is False:
            df = df.query("fold_enrichment>=1").copy()
            logger.info("Found {} GO terms after keeping only positive enrichment".format(len(df)))

        # filter out FDR>0.05
        df = df.query("fdr<=@fdr").copy()
        logger.info("Found {} GO terms after keeping only FDR<{}".format(len(df), fdr))

        return df


    def plot_go_terms(self, ontologies, max_features=50,
        log=False,
        fontsize=8, minimum_genes=0, pvalue=0.05,
        cmap="summer_r",
        sort_by="fold_enrichment", 
        show_pvalues=False,
        include_negative_enrichment=False, 
        fdr_threshold=0.05, compute_levels=True):

        assert sort_by in ['pValue', 'fold_enrichment', 'fdr']

        # FIXME: pvalue and fold_enrichment not sorted in same order
        pylab.clf()

        df = self.get_data(ontologies,
                include_negative_enrichment=include_negative_enrichment,
                fdr=fdr_threshold)

        if len(df) == 0:
            return df


        df = df.query("pValue<=@pvalue")
        logger.info("Filtering out pvalue>{}. Kept {} GO terms".format(pvalue, len(df)))
        df = df.reset_index(drop=True)

        # Select a subset of the data to keep the best max_features in terms of
        # pValue
        subdf = df.query("number_in_list>@minimum_genes").copy()
        logger.info("Filtering out GO terms with less than {} genes: Kept {} GO terms".format(minimum_genes, len(subdf)))

        logger.info("Filtering out the 3 parent terms")
        subdf = subdf.query("id not in @self.ontologies")

        # Keeping only a part of the data, sorting by pValue
        if sort_by == "pValue":
            subdf = subdf.sort_values(by="pValue", ascending=False).iloc[-max_features:]
            df = df.sort_values(by="pValue", ascending=False)
        elif sort_by=="fold_enrichment":
            subdf = subdf.sort_values(by="abs_log2_fold_enrichment", ascending=True).iloc[-max_features:]
            df = df.sort_values(by="abs_log2_fold_enrichment", ascending=False)
        elif sort_by == "fdr":
            subdf = subdf.sort_values(by="fdr", ascending=False).iloc[-max_features:]
            df = df.sort_values(by="fdr", ascending=False)

        subdf = subdf.reset_index(drop=True)

        # We get all levels for each go id.
        # They are stored by MF, CC or BP
        if compute_levels:
            paths = self.get_graph(list(subdf['id'].values))
            levels = []
            keys = list(paths.keys())
            goid_levels = paths[keys[0]]
            if len(keys)>1:
                for k in keys[1:]:
                    goid_levels.update(paths[k])
            levels = [goid_levels[ID] for ID in subdf['id'].values]
            subdf["level"] = levels
        else:
            subdf['level'] = ""
        N = len(subdf)

        size_factor = 12000 / len(subdf)
        max_size = subdf.number_in_list.max()
        min_size = subdf.number_in_list.min()
        sizes = [max(max_size*0.2,x) for x in size_factor * subdf.number_in_list.values/subdf.number_in_list.max()]

        m1  = min(sizes)
        m3  = max(sizes)
        m2  = m1 + (m3-m1)/2

        if log:
            pylab.scatter(pylab.log2(subdf.fold_enrichment), range(len(subdf)),
                c=subdf.fdr,
                s=sizes, 
                cmap=cmap,
                alpha=0.8, 
                ec="k",
                vmin=0, vmax=fdr_threshold, zorder=10)
            #pylab.barh(range(N), pylab.log2(subdf.fold_enrichment), color="r",
            #    label="pvalue>0.05; FDR>0.05")
            #pylab.axvline(1, color="gray", ls="--")
            #pylab.axvline(-1, color="gray", ls="--")
        else:
             pylab.scatter(subdf.fold_enrichment, range(len(subdf)),
                c=subdf.fdr, cmap=cmap, 
                s=sizes, 
                ec="k",alpha=.8, 
                vmin=0, vmax=fdr_threshold, zorder=10)
            #    pylab.barh(range(N), subdf.fold_enrichment, color="r",
            #    label="not significant")
        pylab.grid(zorder=-10)
        ax2 = pylab.colorbar(shrink=0.5); 
        ax2.ax.set_ylabel('FDR')

        

        labels = [x if len(x)<50 else x[0:47]+"..." for x in list(subdf.label)]
        ticks = ["{} ({}) {}".format(ID,level, "; " + label.title()) 
            for level, ID, label in zip(subdf['level'], subdf.id, labels)]

        pylab.yticks(range(N), ticks, fontsize=fontsize,ha='left')

        yax = pylab.gca().get_yaxis()
        try:
            pad =  [x.label.get_window_extent().width for x in yax.majorTicks]
            yax.set_tick_params(pad=max(pad))
        except:
            yax.set_tick_params(pad=60*fontsize*0.7)
        yax.set_tick_params(pad=60*fontsize*0.6)


        fc_max = subdf.fold_enrichment.max(skipna=True)
        fc_min = subdf.fold_enrichment.min(skipna=True)
        # go into log2 space
        fc_max = pylab.log2(fc_max) 
        fc_min = pylab.log2(fc_min) 
        abs_max = max(fc_max, abs(fc_min), 1)

        if log:
            fc_max = abs_max * 1.5
        else:
            fc_max = 2**abs_max * 1.2

        pylab.axvline(0, color="k", lw=2)
        if log:
            pylab.xlabel("Fold Enrichment (log2)")
        else:
            pylab.xlabel("Fold Enrichment")
        if include_negative_enrichment:
            pylab.xlim([-fc_max, fc_max])
        else:
            pylab.xlim([0, fc_max])
        pylab.tight_layout()

        # The pvalue:
        if show_pvalues:
            ax = pylab.gca().twiny()
            ax.set_xlim([0, max(-pylab.log10(subdf.pValue))*1.2])
            ax.set_xlabel("p-values (log10)", fontsize=12)
            ax.plot(-pylab.log10(subdf.pValue), range(len(subdf)), label="pvalue", lw=2, color="k")
            ax.axvline(1.33, lw=1, ls="--", color="grey", label="pvalue=0.05")
            pylab.tight_layout()
            pylab.legend(loc="lower right")
        s1 = pylab.scatter([],[], s=m1, marker='o', color='#555555', ec="k")
        s2 = pylab.scatter([],[], s=m2, marker='o', color='#555555', ec="k")
        s3 = pylab.scatter([],[], s=m3, marker='o', color='#555555', ec="k")

        if len(subdf) <10:
            labelspacing=1.5 * 4
            borderpad=4
            handletextpad=2
        elif len(subdf) <20:
            labelspacing=1.5 * 2
            borderpad=1
            handletextpad=2
        else:
            labelspacing=1.5
            borderpad=2
            handletextpad=2

        if len(subdf)>=3:
            leg = pylab.legend((s1,s2,s3),
                (str(int(min_size)), str(int(min_size + (max_size-min_size)/2)),str(int(max_size))),
                scatterpoints=1,
                loc='lower right',
                ncol=1,
                frameon=True, title="gene-set size",
                labelspacing=labelspacing, borderpad=borderpad, 
                 handletextpad=handletextpad,
                fontsize=8)
        else:
            leg = pylab.legend((s1, ),
                (str(int(min_size)), ),
                scatterpoints=1,
                loc='lower right',
                ncol=1,
                frameon=True, title="gene-set size",
                labelspacing=labelspacing, borderpad=borderpad, 
                 handletextpad=handletextpad,
                fontsize=8)

        frame = leg.get_frame()
        frame.set_facecolor('#b4aeae')
        frame.set_edgecolor('black')
        frame.set_alpha(1)

        self.subdf = subdf
        self.df = df
        return df

    def get_graph(self, go_ids, ontologies=None):
        # Here we filter the data to keep only the relevant go terms as shown in
        # panther pie chart
        import networkx as nx
        gg = nx.DiGraph()

        #assert ontology in ['MF', 'BP', 'CC']
        if ontologies is None:
            ontologies = ['MF', 'BP', 'CC']
        elif isinstance(ontologies, str):
            ontologies = [ontologies]
        ancestors = [self.ancestors[x] for x in ontologies]

        levels = []
        real_ids = []
        obsolets = []
        from easydev import Progress
        pb = Progress(len(go_ids))
        print('Retrieving info for each significant go terms')
        annotations = {}

        for i, go_id in enumerate(go_ids):

            # Some go terms maybe obsolet or renamed. Looking at other functions
            # may not work simply because the ID has changed. 
            info = self.quickgo.get_go_terms(go_id)
            annotations[go_id] = info

            if info[0]['id'] != go_id:
                _id = info[0]['id']
                print('changed {} to {}'.format(go_id, _id))
                annotations[_id] = info
            else:
                _id = go_id
            aspect = info[0]['aspect']
            if info[0]['isObsolete'] is True:
                print("Skipping obsole go terms: {}".format(go_id))
                obsolets.append(go_id)
                continue
            real_ids.append(_id)

            # now figure out the distance to main ancestor
            # we can try several times
            #if _id != self.ancestors[ontology]:
            for ancestor in ancestors:

                edges = self.quickgo.get_go_paths(_id, ancestor)
                if edges == 400:
                    print("Could not retrieve {} to {}".format(_id, ancestor))
                    continue
                if edges["numberOfHits"] == 0:
                    continue
                if len(edges["results"])>=1:
                    for path in edges["results"]:
                        for edge in path:
                            gg.add_edge(edge['child'], edge["parent"]) 
                else:
                    print(_id, edges["results"])
            pb.animate(i+1)

        self.obsolets = obsolets
        self.annotations = annotations
        self.graph = gg
        all_paths = {}
        for ancestor in ancestors:
            if ancestor not in gg:
                continue
            paths = nx.shortest_path_length(gg, target=ancestor)
            for obsolet in obsolets:
                paths[obsolet] = 100
            all_paths[ancestor] = paths

        return all_paths

    def save_chart(self, data, filename="chart.png"):
        # if dataframe, get 'id' column, otherwise expect a list or string of go
        # terms separated by commas
        if isinstance(data, list):
            goids = ",".join(data)
        elif isinstance(data, str):
            goids = data
        elif "id" in data:
            goids = ",".join(list(data['id'].values))

        try:
            goids = [x for x in goids.split(',') if x not in self.obsolets]
        except:
            print("populate obsolets attribute")
        goids = ",".join(goids)
        # remove obsolets

        res = self.quickgo.get_go_chart(goids)
        with open(filename, "wb") as fout:
            fout.write(res.content)

        # 282 ids en entrée, seulement 228 en sortie
        # dans les levels, on a des non molecular function a enlever.


class KeggPathwayEnrichment():
    """DRAFT IN PROGRESS


    Current input is the output of the rnadiff analysis
    ::

        pe = PathwayEnrichment("rnadiff", "eco")
        pe.barplot(pe.enrichment['down'])


        # Save all deregulated pathways found by the enrichment:
        up = pe.save_significant_pathways("up")
        down = pe.save_significant_pathways("down")
        up.to_csv("kegg_pathway_up_regulated.csv")
        down.to_csv("kegg_pathway_down_regulated.csv")

    """
    def __init__(self, folder, organism, alpha=0.05, log2_fc=0):
        print("DRAFT in progress")
        from bioservices import KEGG
        self.kegg = KEGG(cache=True)
        self.kegg.organism = organism

        self.rnadiff = RNADiffResults(folder ,alpha=alpha, log2_fc=log2_fc)
        choices = list(self.rnadiff.gene_lists.keys())

        self.background = len(self.kegg.list(self.kegg.organism).split("\n"))
        logger.info("Set number of genes to {}".format(self.background))

        self._load_pathways()
        self.compute_enrichment()


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
            gene_list = list(self.rnadiff.gene_lists[category])

        enr = gseapy.enrichr(
            gene_list=gene_list,
            gene_sets=self.gene_sets,
            verbose=verbose,
            background=background,
            outdir="test", no_plot=True)

        return enr

    def _get_final_df(self, df, cutoff=0.05, nmax=10):
        # takes the df and populate the name and size of the found pathways
        # we also sort by adjusted p-value
        # we keep adj p-value <=0.05
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

    def save_all_pathways(self): #pragma: no cover
        # This does not do any enrichment. Just save all pathways once for all
        # with useful information
        for ID in self.kegg.pathwayIds:
            self.save_pathway(ID)

    def save_significant_pathways(self, mode, cutoff=0.05, nmax=20,
        background=None):#pragma: no cover  

        if background is None:
            background = self.background

        # select the relevant pathways
        df = self._enrichr(mode, background).results
        df = self._get_final_df(df, cutoff=cutoff, nmax=nmax)
        logger.warning("Found {} pathways to save".format(len(df)))
        if len(df) == nmax:
            logger.warning("Restricted pathways to {}".format(nmax))

        logger.info("saving {} deregulated pathways".format(len(df)))

        summaries = {}
        # save them
        for ID in df['Term']:
            summary = self.save_pathway(ID, filename="{}_{}.png".format(ID, mode))
            summaries[ID] = summary
        return summaries


    def find_pathways_by_gene(self, gene_name, match="exact"):
        """Returns pathways that contain the gene name

        ke.find_pathways_by_gene("ysgA") 
        """

        #First let us find the kegg ID
        genes = self.kegg.list(self.kegg.organism).strip().split("\n")

        keggid = [x.split("\t")[0].strip() for x in genes]
        gene_names = [x.split("\t")[1].split(";")[0].strip() for x in genes]

        self.keggid = keggid
        self.gene_names = gene_names
        candidates = []
        for x,y in zip(keggid, gene_names):

            if match == "exact":
                if gene_name == y:
                    candidates = x.split(":")[1]
                    break
            else:
                if gene_name in y:
                    candidates.append(x)
        if match != "exact":
            candidates = [x.split(":")[1] for x in candidates]
            logger.info("Found {} candidate(s): {}".format(len(candidates), candidates))
        else:
            logger.info("Found {} in {}".format(gene_name, candidates))



        paths = []
        for key in self.pathways.keys():
            if "GENE" in self.pathways[key]:
                if match == "exact":
                    if candidates in self.pathways[key]['GENE'].keys():
                        paths.append(key)
                else:
                    for candidate in candidates:
                        if candidate in self.pathways[key]['GENE'].keys():
                            paths.append(key)
        return list(set(paths))






















