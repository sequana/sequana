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

from pathlib import Path
import os
import json

from sequana.lazy import pandas as pd
from sequana.lazy import pylab
from sequana.lazy import numpy as np
from matplotlib_venn import venn2_unweighted, venn3_unweighted
import gseapy

from sequana.summary import Summary

import colorlog

logger = colorlog.getLogger(__name__)


__all__ = ["PantherEnrichment"]


class PantherEnrichment:
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
        max_entries=2000,
        annot_col="Name",
    ):
        """


        rnadiff if provided, superseeds the input filename. This is useful for
        debugging
        """

        self.gene_lists = gene_lists
        self.enrichment_fdr = enrichment_fdr

        # users can set the fold change threshold in the log2 scale or normal
        # scale.
        assert log2_fc_threshold >= 0, "log2 fc_threshold must be >=0"
        if fc_threshold is not None:
            log2_fc_threshold = pylab.log2(fc_threshold)

        from bioservices import panther, quickgo

        self.panther = panther.Panther(cache=True)
        self.valid_taxons = [
            x["taxon_id"] for x in self.panther.get_supported_genomes()
        ]
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
        self.ontologies = [
            "GO:0003674",
            "GO:0008150",
            "GO:0005575",
            "ANNOT_TYPE_ID_PANTHER_GO_SLIM_MF",
            "ANNOT_TYPE_ID_PANTHER_GO_SLIM_BP",
            "ANNOT_TYPE_ID_PANTHER_GO_SLIM_CC",
            "ANNOT_TYPE_ID_PANTHER_PC",
            "ANNOT_TYPE_ID_PANTHER_PATHWAY",
            "ANNOT_TYPE_ID_REACTOME_PATHWAY",
        ]
        self.MF = "GO:0003674"
        self.CC = "GO:0005575"
        self.BP = "GO:0008150"

        self.ontology_aliases = [
            "MF",
            "BP",
            "CC",
            "SLIM_MF",
            "SLIM_BP",
            "SLIM_CC",
            "PROTEIN",
            "PANTHER_PATHWAY",
            "REACTOME_PATHWAY",
        ]

        # panther accepts onyl ~2-3000 genes at max. Let us restrict the analysis
        # to the first 2000 genes based on their log2 fold change 2000 + and
        # 2000 negatives

        msg = "Ignoring DEGs with adjusted p-value > {} and fold change in [{}, {}]".format(
            padj_threshold, 1 / (2 ** log2_fc_threshold), 2 ** log2_fc_threshold
        )
        logger.info(msg)

        # used in report module
        self.summary["fold_change_range"] = [
            1 / (2 ** log2_fc_threshold),
            2 ** log2_fc_threshold,
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
            raise ValueError(
                f"taxon {taxon} not in pantherDB. please use one of {self.valid_taxons}"
            )
        self.taxon_info = [
            x for x in self.panther.get_supported_genomes() if x["taxon_id"] == taxon
        ]
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
        if isinstance(mygenes, list):
            mygenes = ",".join(mygenes)

        if mygenes.count(",") > 2500:
            logger.warning(
                "Please reduce the list input genes. may fail on pantherb otherwise"
            )
        if len(mygenes) <= 2:
            logger.error(
                f"Less than 2 genes are found for in the gene set: {mygenes}. No enrichment will be computed"
            )
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

        for ontology in ontologies:
            logger.info(f" - Computing enrichment for {ontology}")

            results = self.panther.get_enrichment(
                mygenes,
                taxid,
                get_panther_ont(ontology),
                enrichment_test=enrichment_test,
                correction=correction,
            )
            count = 0
            while count < 2 and results == 404:
                logger.warning("Panther request failed Trying again")
                results = self.panther.get_enrichment(
                    mygenes,
                    taxid,
                    get_panther_ont(ontology),
                    enrichment_test=enrichment_test,
                    correction=correction,
                )
                count += 1

            if results == 404:
                logger.warning(
                    "Invalid output from pantherdb (too many genes ?). skipping {}".format(
                        ontology
                    )
                )
                enrichment[ontology] = None
                continue

            if isinstance(results["result"], dict):  # pragma: no cover
                results["result"] = [results["result"]]
            pvalues = [x["pValue"] for x in results["result"]]
            import statsmodels
            import statsmodels.stats.multitest

            if correction == "FDR":
                fdr = statsmodels.stats.multitest.multipletests(
                    pvalues, method="fdr_bh"
                )[1]
            elif correction.lower() == "bonferroni":
                fdr = statsmodels.stats.multitest.multipletests(
                    pvalues, method="bonferroni"
                )[1]
            for i, pvalue in enumerate(pvalues):
                results["result"][i]["fdr2"] = fdr[i]
                if enrichment_test.lower() == "binomial":
                    results["result"][i]["fdr"] = fdr[i]

            enrichment[ontology] = results
        stats = dict([(k, len(v["result"])) for k, v in enrichment.items()])
        stats["input_genes"] = len(mygenes.split(","))

        try:
            unmapped = enrichment[ontologies[0]]["input_list"]["unmapped_id"]
            stats["unmapped_genes"] = unmapped
            stats["N_unmapped_genes"] = len(unmapped)
        except:
            stats["unmapped_genes"] = []
            stats["N_unmapped_genes"] = 0

        # Here, looking at the FDr, it appears that when using bonferroni,
        # all FDR are set to zeros. Moreover, when using Fisher tests and
        # FDR (supposibly a FDR_BH, the results are noisy as compare to a
        # test from statsmodels. Moreover, when using binomial test, the FDR
        # is not computed... So, we will recompute the FDR ourself
        return enrichment, stats

    def get_functional_classification(
        self, mygenes, taxon
    ):  # pragma: no cover ; too slow
        """Mapping information from pantherDB for the lisf of genes

        We also store uniprot persistent id

        """
        logger.warning("Very slow. Please wait")
        if isinstance(mygenes, list):
            mygenes = ",".join(mygenes)

        res = self.panther.get_mapping(mygenes, taxon)
        res = res["mapped"]
        N = len(res)

        from easydev import Progress

        pb = Progress(N)

        for i, item in enumerate(res):
            accession = item["accession"]
            res[i]["persistent_id"] = self._get_name_given_accession(accession)
            pb.animate(i + 1)
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

    def get_data(
        self, category, ontologies, include_negative_enrichment=True, fdr=0.05
    ):
        """

        From all input GO term that have been found and stored in
        enrichment[ONTOLOGY]['result'], we keep those with fdr<0.05. We also
        exclude UNCLASSIFIED entries. The final dataframe is returned

        ::

            pe.get_data("MF")

        """
        if isinstance(ontologies, str):
            ontologies = [ontologies]
        else:
            assert isinstance(ontologies, list)

        if category not in self.enrichment:
            logger.warning("You must call compute_enrichment_{}".format(category))
            return

        # First, we select the required ontologies and build a common data set
        all_data = []
        for ontology in ontologies:
            data = self.enrichment[category][ontology]["result"]
            if isinstance(data, dict):
                # there was only one hit, we expect:
                data = [data]
            all_data.extend(data)
        data = all_data

        # remove unclassified GO terms
        unclassified = [x for x in data if x["term"]["label"] == "UNCLASSIFIED"]
        logger.info("Found {} unclassified".format(len(unclassified)))
        data = [x for x in data if x["term"]["label"] != "UNCLASSIFIED"]

        df = pd.DataFrame(data)
        if len(df) == 0:
            return df
        else:
            logger.info("Found {} GO terms".format(len(df)))

        df = df.query("number_in_list!=0").copy()
        logger.info(
            "Found {} GO terms with at least 1 gene in reference".format(len(df))
        )

        # extract the ID and label
        df["id"] = [x["id"] for x in df["term"]]
        df["label"] = [x["label"] for x in df["term"]]

        # some extra information for convenience
        df["pct_diff_expr"] = df["number_in_list"] * 100 / df["number_in_reference"]
        df["log2_fold_enrichment"] = pylab.log2(df["fold_enrichment"])
        df["abs_log2_fold_enrichment"] = abs(pylab.log2(df["fold_enrichment"]))
        df["expected"] = [int(x) for x in df.expected]

        # Some user may want to include GO terms with fold enrichment
        # significanyly below 1 or not.
        if include_negative_enrichment is False:
            df = df.query("fold_enrichment>=1").copy()
            logger.info(
                "Found {} GO terms after keeping only positive enrichment".format(
                    len(df)
                )
            )

        # filter out FDR>0.05
        df = df.query("fdr<=@fdr").copy()
        logger.info("Found {} GO terms after keeping only FDR<{}".format(len(df), fdr))

        return df

    def _get_plot_go_terms_data(self, category, ontologies=None,
        max_features=50, minimum_genes=0, pvalue=0.05,
        sort_by="fold_enrichment", fdr_threshold=0.05,
        include_negative_enrichment=False, compute_levels=False):

        if ontologies is None:
            ontologies = ["MF", "BP", "CC"]
        assert sort_by in ["pValue", "fold_enrichment", "fdr"]

        df = self.get_data(
            category,
            ontologies,
            include_negative_enrichment=include_negative_enrichment,
            fdr=fdr_threshold,
        )

        if df is None or len(df) == 0:
            return None, None

        # df stores the entire data set
        # subdf will store the subset (max of n_features, and add dummy values)
        df = df.query("pValue<=@pvalue")
        df = df.reset_index(drop=True)

        # Select a subset of the data to keep the best max_features in terms of
        # pValue
        if minimum_genes != 0:
            subdf = df.query("number_in_list>@minimum_genes").copy()
            logger.info(
                f"Keeping {len(subdf)} GO terms with at least {minimum_genes+1} genes"
            )
        else:
            subdf = df.copy()

        logger.debug("Filtering out the 3 parent terms")
        to_ignore = [self.MF, self.CC, self.BP]
        subdf = subdf.query("id not in @to_ignore")
        df = df.query("id not in @to_ignore")

        if subdf is None or len(subdf) == 0:
            return subdf

        # Keeping only a part of the data, sorting by pValue
        if sort_by == "pValue":
            subdf = subdf.sort_values(by="pValue", ascending=False).iloc[-max_features:]
            df = df.sort_values(by="pValue", ascending=False)
        elif sort_by == "fold_enrichment":
            subdf = subdf.sort_values(
                by="abs_log2_fold_enrichment", ascending=True
            ).iloc[-max_features:]
            df = df.sort_values(by="abs_log2_fold_enrichment", ascending=False)
        elif sort_by == "fdr":
            subdf = subdf.sort_values(by="fdr", ascending=False).iloc[-max_features:]
            df = df.sort_values(by="fdr", ascending=False)

        subdf = subdf.reset_index(drop=True)

        # We get all levels for each go id. They are stored by MF, CC or BP
        subdf["level"] = ""
        if compute_levels:
            paths = self._get_graph(
                list(subdf["id"].values), ontologies=ontologies
            )
            if paths:
                levels = []
                keys = list(paths.keys())

                # FIXME this part is flaky. What would happen if the levels are
                # different if several keys are found ? We use the last one...
                goid_levels = paths[keys[0]]
                if len(keys) > 1:
                    for k in keys[1:]:
                        goid_levels.update(paths[k])

                levels = [goid_levels[ID] for ID in subdf["id"].values]
                subdf["level"] = levels

        return df, subdf

    def plot_go_terms(
        self,
        category,
        ontologies=None,
        max_features=50,
        log=False,
        fontsize=9,
        minimum_genes=0,
        pvalue=0.05,
        cmap="summer_r",
        sort_by="fold_enrichment",
        show_pvalues=False,
        include_negative_enrichment=False,
        fdr_threshold=0.05,
        compute_levels=True,
        progress=True,
    ):

        df, subdf = self._get_plot_go_terms_data(category, ontologies=ontologies,
            max_features=max_features, minimum_genes=minimum_genes, pvalue=pvalue,
            include_negative_enrichment=include_negative_enrichment,
            sort_by=sort_by, fdr_threshold=fdr_threshold,
            compute_levels=compute_levels)

        if df is None or subdf is None:
            return

        # now, for the subdf, which is used to plot the results, we add dummy
        # rows to make the yticks range scale nicer.
        M = 10
        datum = subdf.iloc[-1].copy()
        datum.fdr = 0
        datum.number_in_list = 0
        datum.fold_enrichment = 1
        datum.label = ""
        datum["id"] = ""
        datum["level"] = ""

        while len(subdf) < 10:
            subdf = pd.concat([datum.to_frame().T, subdf], axis=0)

        # here, we try to figure out a proper layout
        N = len(subdf)
        size_factor = 10000 / len(subdf)
        max_size = subdf.number_in_list.max()
        # ignore the dummy values
        min_size = min([x for x in subdf.number_in_list.values if x != 0])
        # here we define a size for each GO entry.
        # For the dummy entries, size is null (int(bool(x))) makes sure
        # it is not shown
        sizes = [
            max(max_size * 0.2, x) * int(bool(x))
            for x in size_factor
            * subdf.number_in_list.values
            / subdf.number_in_list.max()
        ]

        m1 = min([x for x in sizes if x != 0])
        m3 = max(sizes)
        m2 = m1 + (m3 - m1) / 2

        # The plot itself. we stretch when there is lots of features
        if len(subdf) > 25:
            fig = pylab.figure(num=1)
            fig.set_figwidth(10)
            fig.set_figheight(8)
        else:
            fig = pylab.figure(num=1)
            fig.set_figwidth(10)
            fig.set_figheight(6)

        pylab.clf()
        if log:
            pylab.scatter(
                [pylab.log2(x) if x else 0 for x in subdf.fold_enrichment],
                range(len(subdf)),
                c=subdf.fdr,
                s=sizes,
                cmap=cmap,
                alpha=0.8,
                ec="k",
                vmin=0,
                vmax=fdr_threshold,
                zorder=10,
            )
        else:
            pylab.scatter(
                subdf.fold_enrichment,
                range(len(subdf)),
                c=subdf.fdr,
                cmap=cmap,
                s=sizes,
                ec="k",
                alpha=0.8,
                vmin=0,
                vmax=fdr_threshold,
                zorder=10,
            )

        # set color bar height
        pylab.grid(zorder=-10)
        ax2 = pylab.colorbar(shrink=0.5)
        ax2.ax.set_ylabel("FDR")

        # define the labels
        max_label_length = 45
        labels = [
            x if len(x) < max_label_length else x[0 : max_label_length - 3] + "..."
            for x in list(subdf.label)
        ]
        ticks = []
        for level, ID, label in zip(subdf["level"], subdf.id, labels):
            if ID:
                if level:
                    ticks.append(f"{ID} ({level}) ;  {label.title()}")
                else:
                    ticks.append(f"{ID} ; {label.title()}")
            else:
                ticks.append("")

        # Refine the fontsize of ylabel if not many
        if len(subdf) < 10:
            pylab.yticks(range(N), ticks, fontsize=fontsize, ha="left")
        else:
            pylab.yticks(range(N), ticks, fontsize=fontsize, ha="left")

        yax = pylab.gca().get_yaxis()
        try:
            pad = [x.label1.get_window_extent().width for x in yax.majorTicks]
            yax.set_tick_params(pad=max(pad))
        except:
            yax.set_tick_params(pad=60 * fontsize * 0.7)
        yax.set_tick_params(pad=60 * fontsize * 0.6)

        # deal with the x-axis now. what is the range ?
        fc_max = subdf.fold_enrichment.max(skipna=True)
        fc_min = subdf.fold_enrichment.min(skipna=True)
        # go into log2 space
        fc_max = pylab.log2(fc_max)
        fc_min = pylab.log2(fc_min)
        abs_max = max(fc_max, abs(fc_min), 1)

        if log:
            fc_max = abs_max * 1.5
        else:
            fc_max = 2 ** abs_max * 1.2

        pylab.axvline(0, color="k", lw=2)
        if log:
            pylab.xlabel("Fold Enrichment (log2)")
        else:
            pylab.xlabel("Fold Enrichment")

        # dealwith fold change below 0.
        if include_negative_enrichment:
            pylab.xlim([-fc_max, fc_max])
        else:
            pylab.xlim([0, fc_max])
        pylab.tight_layout()

        # The pvalues:
        if show_pvalues:
            ax = pylab.gca().twiny()
            # ax.set_xlim([0, max(-pylab.log10(subdf.pValue))*1.2])
            pvalues = [-pylab.log10(pv) if pv > 0 else 0 for pv in subdf.pValue]

            ax.set_xlim([0, max(pvalues) * 1.2])
            ax.set_xlabel("p-values (log10)", fontsize=12)
            ax.plot(pvalues, range(len(subdf)), label="pvalue", lw=2, color="k")
            ax.axvline(1.33, lw=1, ls="--", color="grey", label="pvalue=0.05")
            pylab.tight_layout()
            pylab.legend(loc="lower right")

        # now, let us add a legend
        s1 = pylab.scatter([], [], s=m1, marker="o", color="#555555", ec="k")
        s2 = pylab.scatter([], [], s=m2, marker="o", color="#555555", ec="k")
        s3 = pylab.scatter([], [], s=m3, marker="o", color="#555555", ec="k")

        if len(subdf) <= 10:
            labelspacing = 1.5 * 2
            borderpad = 1.5
            handletextpad = 2
        elif len(subdf) < 20:
            labelspacing = 1.5 * 2
            borderpad = 1
            handletextpad = 2
        else:
            labelspacing = 1.5
            borderpad = 2
            handletextpad = 2

        # get back the dataframe without the dummies
        subdf = subdf.query("number_in_list>0")
        if len(subdf) >= 3:
            leg = pylab.legend(
                (s1, s2, s3),
                (
                    str(int(min_size)),
                    str(int(min_size + (max_size - min_size) / 2)),
                    str(int(max_size)),
                ),
                scatterpoints=1,
                loc="lower right",
                ncol=1,
                frameon=True,
                title="gene-set size",
                labelspacing=labelspacing,
                borderpad=borderpad,
                handletextpad=handletextpad,
                fontsize=8,
            )
        elif len(subdf) >= 2:
            leg = pylab.legend(
                (s1, s3),
                (str(int(min_size)), str(int(max_size))),
                scatterpoints=1,
                loc="lower right",
                ncol=1,
                frameon=True,
                title="gene-set size",
                labelspacing=labelspacing,
                borderpad=borderpad,
                handletextpad=handletextpad,
                fontsize=8,
            )
        else:
            leg = pylab.legend(
                (s1,),
                (str(int(min_size)),),
                scatterpoints=1,
                loc="lower right",
                ncol=1,
                frameon=True,
                title="gene-set size",
                labelspacing=labelspacing,
                borderpad=borderpad,
                handletextpad=handletextpad,
                fontsize=8,
            )

        frame = leg.get_frame()
        frame.set_facecolor("#b4aeae")
        frame.set_edgecolor("black")
        frame.set_alpha(1)

        return df

    def _get_graph(self, go_ids, ontologies=None):
        # Here we filter the data to keep only the relevant go terms as shown in
        # panther pie chart
        import networkx as nx

        gg = nx.DiGraph()

        if isinstance(ontologies, str):
            ontologies = [ontologies]

        for x in ontologies:
            if "PROTEIN" in x or "PATHWAY" in x:
                return {}

        ancestors = [self._ancestors[x] for x in ontologies]
        levels = []
        renamed_ids = {}
        obsolets = []
        from easydev import Progress

        logger.info(f"Retrieving info for {len(go_ids)} enriched go terms")
        annotations = {}

        for i, go_id in enumerate(go_ids):

            # retrive info about a given GO ID
            info = self.quickgo.get_go_terms(go_id)
            annotations[go_id] = info

            # Some go terms may be renamed.
            if info[0]["id"] != go_id:
                _id = info[0]["id"]
                logger.warning("changed {} to {}".format(go_id, _id))
                annotations[_id] = info
                renamed_ids[go_id] = _id
            else:
                _id = go_id

            # Some go terms may be obsolete
            if info[0]["isObsolete"] is True:
                logger.warning("Skipping obsolete go terms: {}".format(go_id))
                obsolets.append(go_id)
                continue

            # now figure out the distance to main ancestor
            # we can try several times
            # if _id != self.ancestors[ontology]:
            for ancestor in ancestors:
                edges = self.quickgo.get_go_paths(_id, ancestor)
                if edges == 400:
                    logger.warning("Could not retrieve {} to {}".format(_id, ancestor))
                    continue
                if edges["numberOfHits"] == 0:
                    continue
                if len(edges["results"]) >= 1:
                    for path in edges["results"]:
                        for edge in path:
                            gg.add_edge(edge["child"], edge["parent"])
                else:
                    print(_id, edges["results"])

        self.obsolets += obsolets
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

        for key in all_paths.keys():
            for old, new in renamed_ids.items():
                if new in all_paths[key]:
                    all_paths[key][old] = all_paths[key][new]

        return all_paths

    def save_chart(self, data, filename="chart.png"):
        """

        pe = PantherEnrichment("B4052-V1.T1vsT0.complete.xls", fc_threshold=5,
            padj_threshold=0.05)
        df = pe.plot_go_terms("down", log=True, compute_levels=False)
        pe.save_chart(df, "chart.png")

        """
        # if dataframe, get 'id' column, otherwise expect a list or string of go
        # terms separated by commas
        if isinstance(data, list):
            goids = ",".join(data)
        elif isinstance(data, str):
            goids = data
        elif "id" in data:
            goids = ",".join(list(data["id"].values))

        try:
            goids = [x for x in goids.split(",") if x not in self.obsolets]
        except:
            logger.error("Could not save chart")
        goids = ",".join(goids)
        # remove obsolets

        try:
            res = self.quickgo.get_go_chart(goids)

            if res is None:
                raise Exception
            with open(filename, "wb") as fout:
                fout.write(res.content)
        except:
            import shutil

            logger.warning(
                "Could not create the GO chart. Maybe too many go IDs ({})".format(
                    len(goids.split(","))
                )
            )
            from sequana import sequana_data

            no_data = sequana_data("no_data.png")
            shutil.copy(no_data, filename)


"""BOOK keeping
    def to_excel(self):
        with pd.ExcelWriter(out_dir.parent / "enrichment_go.xlsx") as writer:
            df = self.enrichment_go.copy()
            df.reset_index(inplace=True)
            df.to_excel(writer, "go", index=False)
            ws = writer.sheets["go"]
            try:
                ws.autofilter(0, 0, df.shape[0], df.shape[1] - 1)
            except:
                logger.warning("XLS formatting issue.")
"""
