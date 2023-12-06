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
import plotly.express as px

from sequana.lazy import pandas as pd
from sequana.lazy import pylab

logger = colorlog.getLogger(__name__)


__all__ = ["UniProtEnrichment"]


class PlotGOTerms:
    """Used by :class:`sequana.enrichment.panther.PantherEnrichment`
    and :class:`sequana.enrichment.panther.UniprotEnrichment`

    """

    def __init__(self):
        # will be defined by the parent class
        self.enrichment = {}
        self.gene_sets = {}
        self.df_genes = None

    def get_data(self, category, ontologies, include_negative_enrichment=True, fdr=0.05):
        """

        From all input GO term that have been found and stored in
        enrichment[ONTOLOGY]['result'], we keep those with fdr<0.05. We also
        exclude UNCLASSIFIED entries. The final dataframe is returned

        ::

            pe.get_data("up", "MF")

        """
        if isinstance(ontologies, str):
            ontologies = [ontologies]
        else:
            assert isinstance(ontologies, list)

        if category not in self.enrichment:
            logger.warning(f"Category {category} not found. Have you called compute_enrichment ?")
            return

        # First, we select the required ontologies and build a common data set
        all_data = []
        for ontology in ontologies:
            if ontology not in self.enrichment[category]:
                logger.warning(f"Ontology {ontology} not found. Have you called compute_enrichment ?")
                return

            data = self.enrichment[category][ontology]["result"]
            data["ontology"] = ontology
            all_data.append(data)

        df = pd.concat(all_data, axis=0)

        if len(df) == 0:
            return df
        else:
            logger.info("Found {} GO terms".format(len(df)))

        logger.info("Found {} GO terms with at least 1 gene in reference".format(len(df)))

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

    def _get_plot_go_terms_data(
        self,
        category,
        ontologies=None,
        max_features=50,
        minimum_genes=0,
        pvalue=0.05,
        sort_by="fold_enrichment",
        fdr_threshold=0.05,
        include_negative_enrichment=False,
        compute_levels=False,
    ):
        if ontologies is None:
            ontologies = {"MF", "BP", "CC"}
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
        subdf = df.copy()

        logger.debug("Filtering out the 3 parent terms")
        to_ignore = {"GO:0003674", "GO:0008150", "GO:0005575"}
        subdf = subdf.query("id not in @to_ignore")
        df = df.query("id not in @to_ignore")

        if subdf is None or len(subdf) == 0:
            return df, subdf

        # Keeping only a part of the data, sorting by pValue
        if sort_by in ["pValue", "fdr"]:
            subdf = subdf.sort_values(by=sort_by, ascending=False).iloc[-max_features:]
            df = df.sort_values(by=sort_by, ascending=False)
        elif sort_by == "fold_enrichment":
            subdf = subdf.sort_values(by="abs_log2_fold_enrichment", ascending=True).iloc[-max_features:]
            df = df.sort_values(by="abs_log2_fold_enrichment", ascending=False)

        subdf = subdf.reset_index(drop=True)

        # We get all levels for each go id. They are stored by MF, CC or BP
        subdf["level"] = ""
        if compute_levels:
            paths = self._get_graph(list(subdf["id"].values), ontologies=ontologies)
            if paths:
                levels = []
                keys = list(paths.keys())

                # FIXME this part is flaky. What would happen if the levels are
                # different if several keys are found ? We use the last one...
                goid_levels = paths[keys[0]]
                if len(keys) > 1:
                    for k in keys[1:]:
                        goid_levels.update(paths[k])

                # FIXME
                # in rare cases, imd are not found in _get_graph()
                # need to add them back here with a dummy level
                for ID in subdf["id"].values:
                    if ID not in goid_levels:
                        goid_levels[ID] = 10

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
        include_negative_enrichment=False,
        fdr_threshold=0.05,
        compute_levels=True,
        progress=True,
    ):
        df, subdf = self._get_plot_go_terms_data(
            category,
            ontologies=ontologies,
            max_features=max_features,
            minimum_genes=minimum_genes,
            pvalue=pvalue,
            include_negative_enrichment=include_negative_enrichment,
            sort_by=sort_by,
            fdr_threshold=fdr_threshold,
            compute_levels=compute_levels,
        )

        if df is None or subdf is None:
            return

        # Filter out depleted pathway
        if not include_negative_enrichment:
            subdf = subdf.query("fold_enrichment > 1")

        if log:
            subdf["log2_fold_enrichment"] = [pylab.log2(x) if x else 0 for x in subdf.fold_enrichment]
            fig = px.scatter(
                subdf,
                x="log2_fold_enrichment",
                y="label",
                color="fdr",
                size="number_in_list",
                hover_data=["id", "label", "level", "fold_enrichment", "fdr", "number_in_list"],
                color_continuous_scale="Viridis",
                labels={"log2_fold_enrichment": "Fold enrichment (log2)", "label": "GO term", "fdr": "FDR"},
            )
        else:
            fig = px.scatter(
                subdf,
                x="fold_enrichment",
                y="label",
                color="fdr",
                size="number_in_list",
                hover_data=["fold_enrichment", "label", "fdr", "number_in_list", "id", "level"],
                color_continuous_scale="Viridis",
                labels={"log2_fold_enrichment": "Fold enrichment (log2)", "label": "GO term", "fdr": "FDR"},
            )

        # To have all labels displayed
        if len(df) > 20:
            fig.update_layout(height=800)

        return fig

    def save_chart(self, df, filename="chart.png"):
        self.quick_go_graph.save_chart(df, filename)

    def _get_graph(self, df, ontologies):
        return self.quick_go_graph._get_graph(df, ontologies=ontologies)

    def _get_go_description(self, goids):
        return self.quick_go_graph.get_go_description(goids)

    def _get_data(self, category, ontologies, include_negative_enrichment=True, fdr=0.05):
        """

        From all input GO term that have been found and stored in
        enrichment[ONTOLOGY]['result'], we keep those with fdr<0.05. We also
        exclude UNCLASSIFIED entries. The final dataframe is returned

        ::

            pe.get_data("up", "MF")

        """
        if isinstance(ontologies, str):
            ontologies = [ontologies]
        else:
            assert isinstance(ontologies, list)

        if category not in self.enrichment:
            logger.warning(f"Category {category} not found. Have you called compute_enrichment ?")
            return

        # First, we select the required ontologies and build a common data set
        all_data = []
        for ontology in ontologies:
            if ontology not in self.enrichment[category]:
                logger.warning(f"Ontology {ontology} not found. Have you called compute_enrichment ?")
                return

            data = self.enrichment[category][ontology]["result"]
            data["ontology"] = ontology
            all_data.append(data)

        df = pd.concat(all_data, axis=0)

        if len(df) == 0:
            return df
        else:
            logger.info("Found {} GO terms".format(len(df)))

        logger.info("Found {} GO terms with at least 1 gene in reference".format(len(df)))

        return df
