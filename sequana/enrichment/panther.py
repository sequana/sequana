#
#  This file is part of Sequana software
#
#  Copyright (c) 2016-2021 - Sequana Development Team
#
#  Distributed under the terms of the 3-clause BSD license.
#  The full license is in the LICENSE file, distributed with this software.
#
#  website: https://github.com/sequana/sequana
#  documentation: http://sequana.readthedocs.io
#
##############################################################################

import colorlog
from sequana.enrichment.ontology import Ontology
from sequana.enrichment.plot_go_terms import PlotGOTerms
from sequana.enrichment.quickgo import QuickGOGraph
from sequana.lazy import numpy as np
from sequana.lazy import pandas as pd
from sequana.lazy import pylab
from tqdm import tqdm

logger = colorlog.getLogger(__name__)


__all__ = ["PantherEnrichment"]


class PantherEnrichment(Ontology, PlotGOTerms):
    """

    By default, we keep only the genes with a adjusted pvalue <= 0.05. The
    fold_change threshold is on a log2 scale and set to 0 (meaning no
    filtering). Only one value is requested and
    used to filter out positive and negative fold change (on the log2
    scale). In other word, a log2 fold change threshold of 2 means that we
    filter out values between -2 and 2.

    If you prefer to work in natural scale, just set the parameter
    fc_threshold, which overwrite the log2_fc_threshold parameter.

    ::

        pe = PantherEnrichment("input_file.tsv", taxon=10090, log2_fc_threshold=1)
        # compute enrichment for genes down and up
        pe.compute_enrichment()

        # Results for up case is stored in pe.enrichment
        # then, we plot the most mportat go terms
        df_up = pe.plot_go_terms("up")

        df = pe.plot_go_terms("up", pe.MF)
        pe.save_chart(df, "chart_MF_up.png")

        # all 3 main ontology
        df = pe.plot_go_terms("up")
        pe.save_chart(df, "chart_up.png")

    e.stats contains some statistics. One important is the list of unmapped
    genes. The results from the GO enrichment are stored in the attributes
    enrichment. There, we have again adjusted p-value and a fold enrichment,
    which can in turn be filtered or not.

    You can retrieve the cleaned data using the get_data method.

    You can also plot the GO terms that are significantly enriched using::

        e.plot_go_terms(['MF', 'CC', 'BP])

    This function returns the dataframe used during the plotting.

    If you want to look at the up regulated genes only::

        e.compute_enrichment(pe.mygenes["up"], 83333)
        df = e.plot_go_terms(['MF', 'CC', 'BP'],
                log=False, include_negative_enrichment=False,
                fontsize=8, sort_by='fold_enrichment',
                show_pvalues=True, fdr_threshold=0.05)


    The number of genes is limited to about 3100 depending (don't ask me
    why, this seem to be a hard-coded limitation on PantherDB website).
    In such case, you should add a filter e.g on padj or fold change


    """

    def __init__(
        self,
        gene_lists,
        taxon,
        requests_per_sec=10,
        padj_threshold=0.05,
        log2_fc_threshold=0,
        fc_threshold=None,
        enrichment_fdr=0.05,
        annot_col="Name",
    ):
        """


        rnadiff if provided, superseeds the input filename. This is useful for
        debugging
        """
        Ontology.__init__(self)
        PlotGOTerms.__init__(self)

        self.gene_lists = gene_lists
        self.enrichment_fdr = enrichment_fdr

        # users can set the fold change threshold in the log2 scale or normal
        # scale.
        assert log2_fc_threshold >= 0, "log2 fc_threshold must be >=0"
        if fc_threshold is not None:
            log2_fc_threshold = pylab.log2(fc_threshold)

        from bioservices import panther, quickgo

        self.quick_go_graph = QuickGOGraph()

        self.panther = panther.Panther(cache=True)
        self.valid_taxons = [x["taxon_id"] for x in self.panther.get_supported_genomes()]
        self.summary = {}

        self._taxon = None
        self.taxon = taxon

        self.quickgo = quickgo.QuickGO(cache=True)
        self.quickgo.requests_per_sec = requests_per_sec
        self.quickgo.settings.TIMEOUT = 120

        self._ancestors = {
            "MF": "GO:0003674",
            "CC": "GO:0005575",
            "BP": "GO:0008150",
            "SLIM_MF": "GO:0003674",
            "SLIM_CC": "GO:0005575",
            "SLIM_BP": "GO:0008150",
        }
        self.ontologies.extend(
            [
                "ANNOT_TYPE_ID_PANTHER_GO_SLIM_MF",
                "ANNOT_TYPE_ID_PANTHER_GO_SLIM_BP",
                "ANNOT_TYPE_ID_PANTHER_GO_SLIM_CC",
                "ANNOT_TYPE_ID_PANTHER_PC",
                "ANNOT_TYPE_ID_PANTHER_PATHWAY",
                "ANNOT_TYPE_ID_REACTOME_PATHWAY",
            ]
        )

        self.ontology_aliases.extend(
            [
                "SLIM_MF",
                "SLIM_BP",
                "SLIM_CC",
                "PROTEIN",
                "PANTHER_PATHWAY",
                "REACTOME_PATHWAY",
            ]
        )

        # panther accepts onyl ~2-3000 genes at max. Let us restrict the analysis
        # to the first 2000 genes based on their log2 fold change 2000 + and
        # 2000 negatives

        msg = "Ignoring DEGs with adjusted p-value > {} and fold change in [{}, {}]".format(
            padj_threshold, 1 / (2**log2_fc_threshold), 2**log2_fc_threshold
        )
        logger.info(msg)

        # used in report module
        self.summary["fold_change_range"] = [
            1 / (2**log2_fc_threshold),
            2**log2_fc_threshold,
        ]
        self.summary["padj_threshold"] = padj_threshold

        fc_threshold = log2_fc_threshold

        for x in sorted(gene_lists.keys()):

            N = len(gene_lists[x])
            logger.info(f"Starting with {N} genes from category '{x}'")

        self.summary["DGE_after_filtering"] = {k: len(v) for k, v in gene_lists.items()}

        self.enrichment = {}
        self.stats = {}
        self.obsolets = []

    def _set_taxon(self, taxon):
        if taxon not in self.valid_taxons:
            raise ValueError(f"taxon {taxon} not in pantherDB. please use one of {self.valid_taxons}")
        self.taxon_info = [x for x in self.panther.get_supported_genomes() if x["taxon_id"] == taxon]
        self.taxon_info = self.taxon_info[0]
        self._taxon_id = taxon

    def _get_taxon(self):
        return self._taxon_id

    taxon = property(_get_taxon, _set_taxon)

    def compute_enrichment(
        self,
        taxid=None,
        ontologies=None,
        enrichment_test="FISHER",
        correction="FDR",
    ):
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

        for category, gene_list in self.gene_lists.items():
            logger.info(f"Computing enrichment for category {category}")
            self.enrichment[category], self.stats[category] = self._compute_enrichment(
                gene_list,
                taxid=taxid,
                ontologies=ontologies,
                enrichment_test=enrichment_test,
                correction=correction,
            )

    def get_mapping_stats(self):

        results = []
        for k in self.enrichment.keys():
            if self.enrichment[k]:
                ontologies = self.enrichment[k].keys()
                for o in ontologies:
                    M = self.enrichment[k][o]["input_list"]["mapped_count"]
                    U = self.enrichment[k][o]["input_list"]["unmapped_count"]
                    results.append([k, o, M, U])
        df = pd.DataFrame(results)
        df[4] = 100 * df[2] / (df[3] + df[2])
        df[5] = df[2] + df[3]
        df.columns = [
            "category",
            "ontology",
            "mapped",
            "unmapped",
            "mapped_percentage",
            "total",
        ]
        return df

    def _compute_enrichment(
        self,
        mygenes,
        taxid,
        ontologies=None,
        enrichment_test="FISHER",
        correction="FDR",
    ):

        # taxid=83333 # ecoli
        if taxid is None:
            taxid = self.taxon

        # nice bug from panterdb website. if 'php' is a gene name, no process is done
        if 'php' in mygenes:
            mygenes.remove('php')

        if isinstance(mygenes, list):
            mygenes = ",".join(mygenes)

        if mygenes.count(",") > 2500:  # pragma: no cover
            logger.warning("Please reduce the list input genes. may fail on pantherb otherwise")
        if len(mygenes) <= 2:  # pragma: no cover
            logger.error(f"Less than 2 genes are found for in the gene set: {mygenes}. No enrichment will be computed")
            return None, None

        if ontologies is None:
            ontologies = self.ontology_aliases
        else:
            for x in ontologies:
                assert x in self.ontology_aliases, f"{x} not valid ontology"

        # for each ontology, we will store one key/value item
        enrichment = {}

        def get_panther_ont(x):
            index = self.ontology_aliases.index(x)
            return self.ontologies[index]

        unclassified = {}
        for ontology in ontologies:
            unclassified[ontology] = 0
            logger.info(f" - Computing enrichment for {ontology}")

            try:
                del results
            except NameError:
                pass
            results = self.panther.get_enrichment(
                mygenes,
                taxid,
                get_panther_ont(ontology),
                enrichment_test=enrichment_test,
                correction=correction,
            )
            count = 0
            while count < 2 and results == 404:  # pragma: no cover
                logger.warning("Panther request failed Trying again")
                results = self.panther.get_enrichment(
                    mygenes,
                    taxid,
                    get_panther_ont(ontology),
                    enrichment_test=enrichment_test,
                    correction=correction,
                )
                count += 1

            if results == 404:  # pragma: no cover
                logger.warning("Invalid output from pantherdb (too many genes ?). skipping {}".format(ontology))
                enrichment[ontology] = None
                continue

            ##if isinstance(results["result"], dict):  # pragma: no cover
            #    results["result"] = [results["result"]]

            # extract the ID and label of the go term and save
            # as primary keys



            for i, k in enumerate(results["result"]):
                if results["result"][i]["term"]["label"] == "UNCLASSIFIED":
                    unclassified[ontology] += 1
                    # results["result"][i]['label'] = None
                else:
                    results["result"][i]["id"] = k["term"]["id"]
                    results["result"][i]["label"] = k["term"]["label"]

            # convert to dataframe and remove unclassified labels
            df = pd.DataFrame(results["result"])
            df = df[df["label"].isnull() == False]

            # compute pvalues adjusted
            pvalues = df["pValue"].copy()

            import statsmodels
            import statsmodels.stats.multitest

            if correction == "FDR":
                fdr = statsmodels.stats.multitest.multipletests(pvalues, method="fdr_bh")[1]
            elif correction.lower() == "bonferroni":
                fdr = statsmodels.stats.multitest.multipletests(pvalues, method="bonferroni")[1]

            df["fdr2"] = fdr
            if enrichment_test.lower() == "binomial":
                df["fdr"] = fdr

            # store all results
            enrichment[ontology] = results
            enrichment[ontology]["result"] = df

        stats = dict([(k, len(v["result"])) for k, v in enrichment.items()])
        stats["input_genes"] = len(mygenes.split(","))

        try:
            unmapped = enrichment[ontologies[0]]["input_list"]["unmapped_id"]
            stats["unmapped_genes"] = unmapped
            stats["N_unmapped_genes"] = len(unmapped)
        except Exception:
            stats["unmapped_genes"] = []
            stats["N_unmapped_genes"] = 0
        stats["unclassified_GO"] = unclassified

        # Here, looking at the FDr, it appears that when using bonferroni,
        # all FDR are set to zeros. Moreover, when using Fisher tests and
        # FDR (supposibly a FDR_BH, the results are noisy as compare to a
        # test from statsmodels. Moreover, when using binomial test, the FDR
        # is not computed... So, we will recompute the FDR ourself
        return enrichment, stats

    def get_data(self, category, ontologies, include_negative_enrichment=True, fdr=0.05):
        """

        From all input GO term that have been found and stored in
        enrichment[ONTOLOGY]['result'], we keep those with fdr<0.05. We also
        exclude UNCLASSIFIED entries. The final dataframe is returned

        ::

            pe.get_data(""up", MF")

        """

        df = self._get_data(category, ontologies)

        # some extra information for convenience
        df["pct_diff_expr"] = df["number_in_list"] * 100 / df["number_in_reference"]

        # could happen that all fold enrichment are set to 'NaN'
        df = df[df['fold_enrichment'] != 'NaN']
        df["log2_fold_enrichment"] = pylab.log2(df["fold_enrichment"])
        df["abs_log2_fold_enrichment"] = abs(pylab.log2(df["fold_enrichment"]))
        df["expected"] = [int(x) for x in df.expected]

        # Some user may want to include GO terms with fold enrichment
        # significanyly below 1 or not.
        if include_negative_enrichment is False:
            df = df.query("fold_enrichment>=1").copy()
            logger.info("Found {} GO terms after keeping only positive enrichment".format(len(df)))

        # filter out FDR>0.05
        df = df.query("fdr<=@fdr").copy()
        logger.info("Found {} GO terms after keeping only FDR<{}".format(len(df), fdr))

        return df

    def get_functional_classification(self, mygenes, taxon):  # pragma: no cover ; too slow
        """Mapping information from pantherDB for the lisf of genes

        We also store uniprot persistent id

        """
        logger.warning("Very slow. Please wait")
        if isinstance(mygenes, list):
            mygenes = ",".join(mygenes)

        res = self.panther.get_mapping(mygenes, taxon)
        res = res["mapped"]

        for i, item in tqdm(enumerate(res)):
            accession = item["accession"]
            res[i]["persistent_id"] = self._get_name_given_accession(accession)
        return res

    def _get_name_given_accession(self, accession):  # pragma: no cover
        from bioservices import UniProt

        self.uniprot = UniProt(cache=True)
        self.uniprot.requests_per_sec = 10

        acc = [x for x in accession.split("|") if x.startswith("UniProtKB")]
        acc = acc[0].split("=")[1]
        res = self.uniprot.get_df(acc, limit=1)
        name = res["Gene names  (primary )"][0]
        return name

    def plot_piechart(self, df):
        # Here we show the GO terms that have number in list > 0
        # Note, that this is dangerous to look only at this picture without
        # the reference plot, which data is not available thourg the pathner API
        labels = []
        for this in df.query("number_in_list!=0").label.values:
            if len(this) > 50:
                labels.append(this)
            else:
                labels.append(this[0:50] + "...")
        pylab.pie(df.query("number_in_list!=0").number_in_list, labels=labels)
        pylab.tight_layout()

    def save_chart(self, df, filename="chart.png"):
        """

        pe = PantherEnrichment("B4052-V1.T1vsT0.complete.xls", fc_threshold=5,
            padj_threshold=0.05)
        df = pe.plot_go_terms("down", log=True, compute_levels=False)
        pe.save_chart(df, "chart.png")

        """
        self.quick_go_graph.save_chart(df, filename)

    def _get_graph(self, df, ontologies):
        return self.quick_go_graph._get_graph(df, ontologies=ontologies)
