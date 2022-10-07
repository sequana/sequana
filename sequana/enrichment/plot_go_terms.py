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
                for ID in subdf['id'].values:
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
        show_pvalues=False,
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
            for x in size_factor * subdf.number_in_list.values / subdf.number_in_list.max()
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
        labels = [x if len(x) < max_label_length else x[0 : max_label_length - 3] + "..." for x in list(subdf.label)]
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
            fc_max = 2**abs_max * 1.2

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
