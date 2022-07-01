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
import io
from collections import Counter, defaultdict

import colorlog
from bioservices import UniProt
from sequana.enrichment.gsea import GSEA
from sequana.enrichment.ontology import Ontology
from sequana.enrichment.plot_go_terms import PlotGOTerms
from sequana.enrichment.quickgo import QuickGOGraph
from sequana.lazy import pandas as pd
from sequana.lazy import pylab

logger = colorlog.getLogger(__name__)


__all__ = ["UniProtEnrichment"]


class UniprotEnrichment(Ontology, PlotGOTerms):
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

        pe = UniprotEnrichment("input_file.tsv", taxon=10090, log2_fc_threshold=1)
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

    """

    def __init__(
        self,
        gene_lists,
        taxon,
        padj_threshold=0.05,
        log2_fc_threshold=0,
        fc_threshold=None,
        enrichment_fdr=0.05,
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

        self.summary = {}
        self.quick_go_graph = QuickGOGraph()

        self._taxon = None
        self.taxon = taxon

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

        try:
            self.df_genes = self._fill_uniprot_for_taxon()

            self.df_genes = self.df_genes[self.df_genes["Gene names"].isnull() == False]
            self.df_genes.index = self._get_gene_names(self.df_genes)

            BP = self.df_genes["Gene ontology (biological process)"]
            CC = self.df_genes["Gene ontology (cellular component)"]
            MF = self.df_genes["Gene ontology (molecular function)"]

            self.gene_sets = {
                "BP": self.get_go_gene_dict(BP),
                "MF": self.get_go_gene_dict(MF),
                "CC": self.get_go_gene_dict(CC),
            }
        except TypeError:
            logger.critical("Uniprot search failed. You may try again if the unitprot server was done")

    def get_go_gene_dict(self, df):
        results = defaultdict(list)
        for name, value in zip(df.dropna().index, df.dropna().values):
            if value:
                for y in value.split(";"):
                    go = y.split()[-1].replace("[", "").replace("]", "")
                    results[go].append(name)
        return results

    def _get_gene_names(self, conv):
        results = []
        for i, x in enumerate(conv["Gene names"]):
            names = x.split()

            # Some gene names are very species-dependent. 
            # A gene name may also be composed of several valid names

            # for cryptococcus, we keep first name starting with CNAG , otherwise we keep the first word only.
            names = [x for x in names if x.startswith("CNAG_")] + [x for x in names if not x.startswith("CNAG_")]

            # irrespective of the list of names, we only keep the first one
            results.append(names[0])
        return results

    def _fill_uniprot_for_taxon(self):
        columns = ",".join(
            [
                "accession",
                "id",
                "gene_name",
                "gene_primary",
                "gene_synonym",
                "gene_oln",
                "gene_orf",
                "organism_name",
                "organism_id",
                "protein_name",
                "ec",
                "go",
                "go_p",
                "go_f",
                "go_c",
                "go-id",
            ]
        )
        uniprot = UniProt(cache=True, verbose=True)
        df = uniprot.search(f"organism_name:{self.taxon}", frmt="tsv", columns=columns)

        try:
            df = pd.read_csv(io.StringIO(df), sep="\t")
        except TypeError:
            logger.critical(f"UniProt call failed and return code {df}")
            return []

        # populate some information
        orgs = Counter(df["Organism"])
        max_count = 0
        best_guess = None
        for k, v in orgs.items():
            if v > max_count:
                best_guess = k
        self.taxon_name = best_guess

        return df

    def compute_enrichment(
        self,
        ontologies=None,
        background=None,
    ):
        """
        :param enrichment_test: Fisher or Binomial
        :param correction: FDR or Bonferonni

        The field **number_in_reference** indicates from the reference, the number
        of genes that have a given ontology term. For instance, 998 genes have
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
            self.enrichment[category], self.stats[category] = self._compute_enrichment(
                gene_list, ontologies=ontologies, background=background
            )

    def _compute_enrichment(
        self,
        mygenes,
        ontologies=None,
        background=None,
    ):

        if background is None:
            background = len(self.df_genes)

        if ontologies is None:
            ontologies = self.ontology_aliases
        else:
            for x in ontologies:
                assert x in self.ontology_aliases, f"{x} not valid ontology"

        # for each ontology, we will store one key/value item
        enrichment = {}
        stats = {}

        for ontology in ontologies:

            gs = GSEA(self.gene_sets[ontology])
            enr = gs.compute_enrichment(mygenes, verbose=False, background=background)
            result = enr.results[enr.results["Adjusted P-value"] < 0.05].copy()
            description = self._get_go_description(result["Term"].values)
            result["description"] = description

            enrichment[ontology] = {"result": result}

        return enrichment, stats

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
        self.quick_go_graph.save_chart(df, filename)

    def _get_graph(self, df, ontologies):
        return self.quick_go_graph._get_graph(df, ontologies=ontologies)

    def _get_go_description(self, goids):
        return self.quick_go_graph.get_go_description(goids)

    def get_data(self, category, ontologies, include_negative_enrichment=True, fdr=0.05):
        """

        From all input GO term that have been found and stored in
        enrichment[ONTOLOGY]['result'], we keep those with fdr<0.05. We also
        exclude UNCLASSIFIED entries. The final dataframe is returned

        ::

            pe.get_data("up", "MF")

        """
        df = self._get_data(category, ontologies)

        df.rename({"P-value": "pValue"}, axis=1, inplace=True)
        df.rename({"Odds Ratio": "fold_enrichment"}, axis=1, inplace=True)
        df.rename({"Adjusted P-value": "fdr"}, axis=1, inplace=True)
        df["id"] = df["Term"]

        # extract the ID and label
        df["label"] = df["description"]

        """The field **number_in_reference** indicates from the reference, the number
        of genes that have a given ontology term. For instance, 998 genes have
        the term. This is stored in **number_in_reference**. If the reference
        contains 4391 genes, and you provided 49
        genes , the **expected** number of genes that have this ontology term is
        49*998/4391 that is 11.1369, which is stored in **"expected**.
        """

        nl = []
        terms = df["Term"]
        ontologies = df["ontology"]
        for term, ontology in zip(terms, ontologies):
            nl.append(len(self.gene_sets[ontology][term]))
        df["number_in_reference"] = nl
        df["number_in_list"] = len(df["Genes"])
        df["total_genes"] = len(self.df_genes)

        logger.warning("fold enrichment currently set to log10(adjusted pvalue)")
        df["fold_enrichment"] = -pylab.log10(df["pValue"])

        df["log2_fold_enrichment"] = pylab.log2(df["fold_enrichment"])
        df["abs_log2_fold_enrichment"] = abs(pylab.log2(df["fold_enrichment"]))

        # filter out FDR>0.05
        df = df.query("fdr<=@fdr").copy()
        logger.info("Found {} GO terms after keeping only FDR<{}".format(len(df), fdr))

        return df
