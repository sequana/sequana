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
import glob
import requests

from sequana.lazy import pandas as pd
from sequana.lazy import pylab
from sequana.lazy import numpy as np
from matplotlib_venn import venn2_unweighted, venn3_unweighted
import colormap

import gseapy

from sequana.summary import Summary

import colorlog

logger = colorlog.getLogger(__name__)


__all__ = ["KEGGPathwayEnrichment"]


class KEGGPathwayEnrichment:
    """KEGG Pathways enrichment from DGE results

    DGE = Differentially Gene Expression

    Current input is the output of the RNADiff analysis. This is a file
    than can be read by RNADiffResults

    When performing a GDE analysis, feature counts are computed using an input
    GFF. Depending on your parameters the gene names may be saved as ensembl
    identifiers or gene names. If you have gene names understood by KEGG, you
    simply need to use this code::

        ke = KEGGPathwayEnrichment("rnadiff", "eco") #"eco" for E. coli here

    this calls ke.compute_enrichment() that stores the up, down and all results
    in the attribute :attr:`enrichment` as a dictionary.

    You can now plot the results::

        ke.barplot('down')

    and save enriched pathways as follows::

        up = ke.save_significant_pathways("up")
        down = ke.save_significant_pathways("down")
        up.to_csv("kegg_pathway_up_regulated.csv")
        down.to_csv("kegg_pathway_down_regulated.csv")

    This class works like a charm for ecoli with GFF that uses gene names.

    For mus musculus, organism is **mmu** (not **mus*). you will need to have
    a mapping of the Ensembl ID into KEGG IDs (actually gene name).

    You can perform the conversion using BioServices/BioMart. We have
    implemented a simple function inside Sequana::

        from sequana import Mart
        conv = Mart("mmusculus_gene_ensembl")
        df = conf.query()
        conf.save(df)

    You can then import the dataframe, as follows using the mapper argument::

        import pandas as pd
        df = pd.read_csv("biomart.csv")
        df = df.rename({"external_gene_name":"name", "ensembl_gene_id": "ensembl"},
                axis=1)
        df = df.set_index("ensembl", inplace=True)

        KEGGPathwayEnrichment("path_to_rnadiff", "mmu", mapper=df)

    More generally, when starting KEGGPathwayEnrichment, we read all pathways.

    This may change with time. So, you can save the pathways::

        ke.export_pathways_to_json()

    And read them back::

        ke = KEGGPathwayEnrichment("path_to_rnadiff", "mmu", mapper=df,
            preload_directory="kegg_pathways/mmu")

        df = ke.scatterplot('down')
        tight_layout()
        savefig("B4052_T1vsT0_KE_scatterplot_down.png")
        df = ke.scatterplot('up')
        savefig("B4052_T1vsT0_KE_scatterplot_up.png")

    """

    def __init__(
        self,
        gene_lists,
        organism,
        progress=True,
        mapper=None,
        background=None,
        preload_directory=None,
        convert_input_gene_to_upper_case=False,
    ):
        """


        In some cases, the input identifiers are converted into names thanks to
        the input mapper (csv file). Yet, if the external name are from one
        species and you use another species in kegg, the kegg names may be upper case
        while your species' name are in lower case. In such situations, you may
        set input identifiers are upper case setting the
        convert_input_gene_to_upper_case parameter to True
        """
        self.convert_input_gene_to_upper_case = convert_input_gene_to_upper_case

        from bioservices import KEGG

        self.kegg = KEGG(cache=True)
        self.kegg.organism = organism
        self.summary = Summary("KEGGPathwayEnrichment")
        self.summary.add_params(
            {
                "organism": organism,
                "mapper": (True if mapper is not None else False),
                "background": background,
            }
        )

        self.gene_lists = gene_lists

        if background:
            self.background = background
        else:
            self.background = len(self.kegg.list(self.kegg.organism).split("\n"))
        logger.info("Set number of genes to {}".format(self.background))

        self._load_pathways(progress=progress, preload_directory=preload_directory)

        if isinstance(mapper, str):

            df = pd.read_csv(mapper)
            df = df.rename(
                {"external_gene_name": "name", "ensembl_gene_id": "ensembl"}, axis=1
            )
            df.set_index("ensembl", inplace=True)
            self.mapper = df

        else:  # the dataframe should already contain the correct columns and index
            self.mapper = mapper

        try:
            self.compute_enrichment()
        except Exception as err:
            print(err)
            logger.critical("An error occured while computing enrichments. ")

    def _check_category(self, cat):
        if cat not in self.gene_lists.keys():
            raise ValueError(
                f"category must be set to one of {self.gene_lists.keys()}. You provided {cat}"
            )

    def _load_pathways(self, progress=True, preload_directory=None):
        # This is just loading all pathways once for all
        self.pathways = {}
        if preload_directory:
            # preload is a directory with all pathways in it

            pathways = glob.glob(preload_directory + "/*json")
            for i, name in enumerate(pathways):
                key = name.strip(".json").split("/")[-1]
                with open(name, "r") as fin:
                    data = json.load(fin)
                    self.pathways[key] = data
        else:
            logger.info("loading all pathways from KEGG. may take time the first time")
            from easydev import Progress

            pb = Progress(len(self.kegg.pathwayIds))
            for i, ID in enumerate(self.kegg.pathwayIds):
                self.pathways[ID.replace("path:", "")] = self.kegg.parse(
                    self.kegg.get(ID)
                )
                if progress:
                    pb.animate(i + 1)

        # Some cleanup. Note that if we read the json file, this is different
        # since already cleanup but this code does no harm
        for ID in self.pathways.keys():
            name = self.pathways[ID]["NAME"]
            if isinstance(name, list):
                name = name[0]
            self.pathways[ID]["NAME"] = name.split(" - ", 1)[0]

        # save gene sets
        self.gene_sets = {}
        for ID in self.pathways.keys():
            res = self.pathways[ID]
            if "GENE" in res.keys():
                results = []
                # some pathways reports genes as a dictionary id:'gene name; description' ('.eg. eco')
                # others reports genes as a dictionary id:'description'
                for geneID, description in res["GENE"].items():
                    if ";" in description:
                        name = description.split(";")[0]
                    else:
                        name = geneID
                    results.append(name)

                self.gene_sets[ID] = results
            else:
                logger.debug("SKIPPED (no genes) {}: {}".format(ID, res["NAME"]))

        # save all pathways info
        self.df_pathways = pd.DataFrame(self.pathways).T
        del self.df_pathways["ENTRY"]
        del self.df_pathways["REFERENCE"]
        go = [
            x["GO"] if isinstance(x, dict) and "GO" in x.keys() else None
            for x in self.df_pathways.DBLINKS
        ]
        self.df_pathways["GO"] = go
        del self.df_pathways["DBLINKS"]

    def plot_genesets_hist(self, bins=20):
        N = len(self.gene_sets.keys())
        pylab.clf()
        pylab.hist([len(v) for k, v in self.gene_sets.items()], bins=bins, lw=1, ec="k")
        pylab.title("{} gene sets".format(N))
        pylab.xlabel("Gene set sizes")
        pylab.grid(True)
        a, b = pylab.xlim()
        pylab.xlim([0, b])

    def compute_enrichment(self, background=None):
        if background is None:
            background = self.background

        self.summary.data["missing_genes"] = {}
        self.summary.data["input_gene_list"] = {}

        self.enrichment = {
            category: self._enrichr(category, background=background)
            for category in self.gene_lists.keys()
        }

    def _enrichr(self, category, background=None, verbose=True):

        if background is None:
            background = self.background

        if isinstance(category, list):
            gene_list = category
        else:
            self._check_category(category)
            gene_list = self.gene_lists[category]

        logger.info("Input gene list of {} ids".format(len(gene_list)))
        self.summary.data["input_gene_list"][category] = len(gene_list)

        if self.mapper is not None:
            missing = [x for x in gene_list if x not in self.mapper.index]
            logger.info("Missing genes from mapper dataframe: {}".format(len(missing)))
            self.summary.data["missing_genes"][category] = ",".join(missing)
            gene_list = [x for x in gene_list if x in self.mapper.index]
            identifiers = self.mapper.loc[gene_list]["name"].drop_duplicates().values

            if self.convert_input_gene_to_upper_case:
                identifiers = [x.upper() for x in identifiers if isinstance(x, str)]
            logger.info("Mapped gene list of {} ids".format(len(identifiers)))
            gene_list = list(identifiers)

        enr = gseapy.enrichr(
            gene_list=gene_list,
            gene_sets=self.gene_sets,
            verbose=verbose,
            background=background,
            outdir="test",
            no_plot=True,
        )

        enr.results["Genes"] = [
            ";".join(sorted(x.split(";"))) for x in enr.results["Genes"].values
        ]
        return enr

    def _get_final_df(self, df, cutoff=0.05, nmax=10):
        # takes the df and populate the name and size of the found pathways
        # we also sort by adjusted p-value
        # we keep adj p-value <=0.05
        if len(df) == 0:
            return df

        df = df.copy()
        df["name"] = [self.pathways[x]["NAME"] for x in df.Term]
        df["size"] = [len(x.split(";")) for x in df.Genes]
        df = df.sort_values("Adjusted P-value")
        df.reset_index(drop=True, inplace=True)
        df = df[df["Adjusted P-value"] <= cutoff]

        if len(df) < nmax:
            nmax = len(df)
        df = df.iloc[0:nmax]
        df = df.sort_values("Adjusted P-value", ascending=False)
        df = df.rename({"Term": "pathway_id"}, axis=1)
        df = df[df.columns]
        return df

    def barplot(self, category, cutoff=0.05, nmax=10):
        self._check_category(category)
        df = self._get_final_df(
            self.enrichment[category].results, cutoff=cutoff, nmax=nmax
        )
        if len(df) == 0:
            return df

        pylab.clf()
        pylab.barh(range(len(df)), -pylab.log10(df["Adjusted P-value"]))
        pylab.yticks(range(len(df)), df.name)
        pylab.axvline(1.3, lw=2, ls="--", color="r")
        pylab.grid(True)
        pylab.xlabel("Adjusted p-value (log10)")
        pylab.ylabel("Gene sets")
        a, b = pylab.xlim()
        pylab.xlim([0, b])
        pylab.tight_layout()
        return df

    def scatterplot(self, category, cutoff=0.05, nmax=15, gene_set_size=[]):
        self._check_category(category)
        df = self._get_final_df(
            self.enrichment[category].results, cutoff=cutoff, nmax=nmax
        )
        if len(df) == 0:
            return df

        pylab.clf()
        pylab.scatter(
            -pylab.log10(df["Adjusted P-value"]),
            range(len(df)),
            s=10 * df["size"],
            c=-pylab.log10(df["Adjusted P-value"]),
        )

        pylab.xlabel("Odd ratio")
        pylab.ylabel("Gene sets")

        names = [x[0:40] + "..." if len(x) > 40 else x for x in df.name]

        pylab.yticks(range(len(df)), names)

        a, b = pylab.xlim()
        pylab.xlim([0, b * 1.1 + 1])

        a, b = pylab.ylim()
        pylab.ylim([a - 0.5, b + 0.5])

        pylab.grid(True)
        ax = pylab.gca()

        M = max(df["size"])
        if M > 100:
            l1, l2, l3 = "10", "100", str(M)
        else:
            l1, l2, l3 = str(round(M / 3)), str(round(M * 2 / 3)), str(M)

        handles = [
            pylab.Line2D([0], [0], marker="o", markersize=5, label=l1, ls=""),
            pylab.Line2D([0], [0], marker="o", markersize=10, label=l2, ls=""),
            pylab.Line2D([0], [0], marker="o", markersize=15, label=l3, ls=""),
        ]
        ax.legend(handles=handles, loc="upper left", title="gene-set size")

        pylab.axvline(1.3, lw=2, ls="--", color="r")
        ax = pylab.colorbar(pylab.gci())
        pylab.tight_layout()
        return df

    def _get_summary_pathway(self, pathway_ID, df, l2fc=0, padj=0.05):
        genes = self.df_pathways.loc[pathway_ID]["GENE"]
        df_down = df.query("padj<=@padj and log2FoldChange<@l2fc").copy()
        df_up = df.query("padj<=@padj and log2FoldChange>=@l2fc").copy()

        if "Name" in df_down.columns:
            pass
        elif "gene_name" in df_down.columns:
            df_down["Name"] = df_down["gene_name"]
        elif "ID" in df_down.columns:
            df_down["Name"] = df_down["ID"]
        else:
            raise ValueError("Expected to find a column Name, gene_name or ID")

        if "Name" in df_up.columns:
            pass
        elif "gene_name" in df_up.columns:
            df_up["Name"] = df_up["gene_name"]
        elif "ID" in df_up.columns:
            df_up["Name"] = df_up["ID"]
        else:
            raise ValueError("Expected to find a column Name, gene_name or ID")

        logger.info("{}".format(pathway_ID))
        # logger.info("Total down-regulated: {}".format(len(df_down)))
        # logger.info("Total up-regulated: {}".format(len(df_up)))

        mapper = {}
        for k, v in genes.items():
            mapper[v.split(";")[0]] = k
        self.genes = genes
        self.df_down = df_down
        self.df_up = df_up
        summary_names = []
        summary_keggids = []
        summary_types = []
        summary_pvalues = []
        summary_fcs = []

        if self.mapper is not None:
            if "Name" not in df_down.columns:
                df_down["Name"] = df_down["ID"]
                Names = []
                for index in df_down.index:
                    Names.append(self.mapper.loc[index]["name"][0])
                df_down["Name"] = Names

            if "Name" not in df_up.columns:
                df_up["Name"] = df_up["ID"]
                Names = []
                for index in df_up.index:
                    Names.append(self.mapper.loc[index]["name"][0])
                df_up["Name"] = Names

            #
            identifiers = []
            new_mapper = {}
            for name, kegg_id in mapper.items():
                try:
                    identifier = (
                        self.mapper.query("name == @name")["name"]
                        .drop_duplicates()
                        .index[0]
                    )
                    identifiers.append(identifier)
                    new_mapper[identifier] = kegg_id
                except:
                    logger.warning(
                        "Skipped {}(kegg ID {}). could not find mapping".format(
                            name, kegg_id
                        )
                    )
            mapper = new_mapper

        for name, kegg_id in mapper.items():
            summary_names.append(name)
            summary_keggids.append(kegg_id)

            if name.lower() in {x.lower() for x in df_down.Name.fillna("undefined")}:
                padj = -pylab.log10(df_down.query("Name==@name").padj.values[0])
                fc = df_down.query("Name==@name").log2FoldChange.values[0]
                summary_fcs.append(fc)
                summary_pvalues.append(padj)
                summary_types.append("-")
            elif name.lower() in {x.lower() for x in df_up.Name.fillna("undefined")}:
                padj = -pylab.log10(df_up.query("Name==@name").padj.values[0])
                summary_pvalues.append(padj)
                fc = df_up.query("Name==@name").log2FoldChange.values[0]
                summary_fcs.append(fc)
                summary_types.append("+")
            else:
                summary_pvalues.append(None)
                summary_fcs.append(None)
                summary_types.append("=")

        summary = pd.DataFrame(
            {
                "type": summary_types,
                "name": summary_names,
                "padj": summary_pvalues,
                "fc": summary_fcs,
                "keggid": summary_keggids,
            }
        )
        summary["description"] = [
            self.pathways[pathway_ID]["GENE"][x] for x in summary.keggid
        ]
        return summary

    def _get_colors(self, summary):
        colors = {}
        # replace iterrows by faster method
        for _, row in summary.iterrows():
            row = row.fillna(0)
            l2fc = row["fc"]
            kegg_id = row["keggid"]

            # Orange to Blue using _r means Blue to Orange where Blue is cold
            # (negative fold change)

            cmap = colormap.get_cmap.cmap_builder("PuOr_r")

            if l2fc <= -4:
                l2fc = -4
            elif l2fc >= 4:
                l2fc = 4
            # to get a value in the range 0,1 expectd by the colormap
            l2fc = (l2fc + 4) / 8

            color = colormap.rgb2hex(*cmap(l2fc)[0:3], normalised=True)

            if l2fc < 0.2:
                colors[kegg_id] = f"{color},grey"
            else:
                colors[kegg_id] = f"{color},black"
        return colors

    def save_pathway(self, pathway_ID, df, scale=None, show=False, filename=None):

        summary = self._get_summary_pathway(pathway_ID, df)
        colors = self._get_colors(summary)

        logger.info("pathway {} total genes: {}".format(pathway_ID, len(summary)))
        count_up = len(summary.query("type == '+'"))
        count_down = len(summary.query("type == '-'"))
        logger.info("this pathway down-regulared genes: {}".format(count_down))
        logger.info("this pathway up-regulated genes: {}".format(count_up))

        url = "https://www.kegg.jp/kegg-bin/show_pathway"
        # dcolor = "white"  --> does not work with the post requests unlike get
        # requests
        params = {
            "map": pathway_ID,
            "multi_query": "\r\n".join(
                ["{} {}".format(k, v) for k, v in colors.items()]
            ),
        }

        self.params = params

        html_page = requests.post(url, data=params)

        self.tmp = html_page
        html_page = html_page.content.decode()

        links_to_png = [
            x for x in html_page.split() if "png" in x and x.startswith("src")
        ]
        link_to_png = links_to_png[0].replace("src=", "").replace('"', "")
        r = requests.get("https://www.kegg.jp/{}".format(link_to_png))

        if filename is None:
            filename = "{}.png".format(pathway_ID)

        with open(filename, "wb") as fout:
            fout.write(r.content)

        return summary

    def save_significant_pathways(
        self, category, cutoff=0.05, nmax=20, background=None, tag="", outdir="."
    ):  # pragma: no cover
        """category should be up, down or all"""

        if background is None:
            background = self.background

        # select the relevant pathways
        df = self._enrichr(category, background).results
        df = self._get_final_df(df, cutoff=cutoff, nmax=nmax)
        logger.warning("Found {} pathways to save".format(len(df)))
        if len(df) == nmax:
            logger.warning("Restricted pathways to {}".format(nmax))

        logger.info("saving {} deregulated pathways".format(len(df)))

        summaries = {}

        for ID in df["pathway_id"]:
            summary = self.save_pathway(
                ID, filename=(Path(outdir) / f"{ID}_{category}.png")
            )
            summaries[ID] = summary

        return summaries

    def find_pathways_by_gene(self, gene_name, match="exact"):
        """Returns pathways that contain the gene name

        ke.find_pathways_by_gene("ysgA")
        """

        # First let us find the kegg ID
        genes = self.kegg.list(self.kegg.organism).strip().split("\n")

        keggid = [x.split("\t")[0].strip() for x in genes]
        gene_names = [x.split("\t")[1].split(";")[0].strip() for x in genes]

        self.keggid = keggid
        self.gene_names = gene_names
        candidates = []
        for x, y in zip(keggid, gene_names):

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
                    if candidates in self.pathways[key]["GENE"].keys():
                        paths.append(key)
                else:
                    for candidate in candidates:
                        if candicodate in self.pathways[key]["GENE"].keys():
                            paths.append(key)
        return list(set(paths))

    def save_project(self, tag, outdir="."):
        """Save tables and visualisations of the complete enrichmment analysis."""
        outdir = Path(outdir)
        outdir.mkdir(parents=True, exist_ok=True)

        for category in ["up", "down", "all"]:
            common_out = Path(f"{tag}_kegg_gsea_{category}_degs")

            results = self.enrichment[category].results

            if not results.empty:

                # FIXME: For now fixing a nmax to 10000 to be sure to have all results
                # (this could be improved by having self.complete_df and self.filtered_df attributes)
                self._get_final_df(results, cutoff=1, nmax=10000).to_csv(
                    outdir / (common_out.name + ".csv")
                )
                self._get_final_df(results, cutoff=0.05, nmax=10000).to_csv(
                    outdir / (common_out.name + "_significant.csv")
                )

                self.barplot(category)
                pylab.savefig(outdir / (common_out.name + "_barplot.png"), dpi=200)
                self.scatterplot(category)
                pylab.savefig(outdir / (common_out.name + "_scatterplot.png"), dpi=200)

            # In case of no enrichment results, create empty files stating so
            else:
                (outdir / (common_out.name + "_NO_RESULTS")).touch()

    def export_pathways_to_json(self, outdir="kegg_pathways"):
        # This is useful to keep an exact track of the pathways that were used.
        # They can be loaded back. If so, we use kegg service only in
        # :meth:`find_pathways_by_gene` method and

        outdir = outdir + "/" + self.kegg.organism
        os.makedirs(outdir, exist_ok=True)

        for key, data in self.pathways.items():
            with open(f"{outdir}/{key}.json", "w") as fout:
                json.dump(data, fout)


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
