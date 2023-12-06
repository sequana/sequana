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

import glob
import json
import os
import shutil
import time
import urllib
from os.path import expanduser
from pathlib import Path

import colorlog
import colormap
import plotly.express as px
import requests
from bioservices import KEGG
from tqdm import tqdm

from sequana.enrichment.gsea import GSEA
from sequana.lazy import numpy as np
from sequana.lazy import pandas as pd
from sequana.lazy import pylab
from sequana.summary import Summary

logger = colorlog.getLogger(__name__)


__all__ = ["KEGGPathwayEnrichment"]


class KEGGPathwayEnrichment:
    """KEGG Pathways enrichment from DGE results

    DGE = Differentially Gene Expression

    Current input is the output of the RNADiff analysis. This is a file
    than can be read by RNADiffResults

    This class takes as input a dictionary with 3 lists of genes: 'up',
    'down' and 'all'.

    For example, using the output of the RNA-seq DGE, use::

        from sequana import RNADiffResults
        r = RNADiffResults("rnadiff.csv")
        r._log2_fc = 1
        gl = r.get_gene_lists(annot_col="gene_name")
        gl['KO_vs_cont"]

        ke = KEGGPathwayEnrichment(gl, "mmu") # mmu for mouse

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

        ke = KEGGPathwayEnrichment("path_to_rnadiff", "mmu", mapper=df)

        ke.scatterplot('down')
        savefig("B4052_T1vsT0_KE_scatterplot_down.png")
        ke.scatterplot('up')
        savefig("B4052_T1vsT0_KE_scatterplot_up.png")

    Pathways are loaded from KEGG which may take some time. For development or production, you may want
    to save the KEGG pathways in a given place. Note however, that the pathways that are enriched will
    still be downloaded for the final annotation and image creation

    To save the pathways locally and load them later do as follows (here for human)::

        ke = KEGGPathwayEnrichment({}, organism="hsa")
        ke.save_pathways("all_pathways/")

        # load them back next time
        ke = KEGGPathwayEnrichment({}, organism="hsa", preload_directory="all_pathways")

    """

    def __init__(
        self,
        gene_lists,
        organism,
        progress=True,
        mapper=None,
        background=None,
        padj_cutoff=0.05,
        preload_directory=None,
        convert_input_gene_to_upper_case=False,
        color_node_with_annotation="Name",
        used_genes=None,
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
        self.color_node_with_annotation = color_node_with_annotation

        self.kegg = KEGG(cache=True)
        self.kegg.organism = organism
        self.summary = Summary("KEGGPathwayEnrichment")
        self.gene_lists = gene_lists

        # What is the intersection between the gene names found in KEGG and
        # the used names used within the RNAdiff analysis (and therefore the
        # gene names to be found in gene_lists)
        # First, we get the list of KEGG names and print number of genes
        kegg_gene_names = self.kegg.list(self.kegg.organism).strip()
        N_kegg_genes = len(kegg_gene_names.split("\n"))
        logger.info(f"Number of KEGG genes for {organism}: {N_kegg_genes}")
        # then, we extract the fourth columns that contains gene names and synonyms
        kegg_gene_names = [x.split("\t")[3] for x in kegg_gene_names.split("\n") if len(x.split("\t")) == 4]
        kegg_gene_names = [y.strip() for x in kegg_gene_names for y in x.split(";")]
        kegg_gene_names = [y.strip() for x in kegg_gene_names for y in x.split(",")]

        def _get_intersection_kegg_genes(genes):
            if len(genes) == 0:
                return "NA"
            else:
                return len(set(kegg_gene_names).intersection(set(genes))) / len(genes) * 100

        dge_genes = set([x for k, v in gene_lists.items() for x in v])

        mapped = [
            _get_intersection_kegg_genes(gene_lists["down"]),
            _get_intersection_kegg_genes(gene_lists["up"]),
            _get_intersection_kegg_genes(gene_lists["all"]),
        ]

        lengths = [len(gene_lists["down"]), len(gene_lists["up"]), len(gene_lists["all"])]
        category = ["down", "up", "all"]

        if used_genes is not None:
            lengths.append(len(used_genes))
            category.append("all genes")
            intersection = _get_intersection_kegg_genes(used_genes)
            logger.info(f"Percentage of genes used in RNAdiff and found in KEGG: {intersection}")
            mapped.append(intersection)
            # if no background provided, let us use the number of found genes
            if background is None:
                background = len(used_genes)

        self.overlap_stats = pd.DataFrame({"category": category, "mapped_percentage": mapped, "N": lengths})

        # Define the background
        if background:
            self.background = background
            logger.info(f"User-defined background: {self.background}")
        else:
            self.background = N_kegg_genes
            logger.info(f"Background defined as number of KEGG genes: {self.background} ")

        self.summary.add_params(
            {
                "organism": organism,
                "mapper": (True if mapper is not None else False),
                "background": self.background,
            }
        )

        self.padj_cutoff = padj_cutoff

        self._load_pathways(progress=progress, preload_directory=preload_directory)

        if isinstance(mapper, str):
            df = pd.read_csv(mapper)
            df = df.rename({"external_gene_name": "name", "ensembl_gene_id": "ensembl"}, axis=1)
            df.set_index("ensembl", inplace=True)
            self.mapper = df

        else:  # the dataframe should already contain the correct columns and index
            self.mapper = mapper

        try:
            self.compute_enrichment()
        except Exception as err:  # pragma: no cover
            logger.critical(err)
            logger.critical("An error occured while computing enrichments. ")

    def _check_category(self, cat):
        if cat not in self.gene_lists.keys():
            raise ValueError(f"category must be set to one of {self.gene_lists.keys()}. You provided {cat}")

    def save_pathways(self, out_directory):
        outdir = Path(out_directory)
        outdir.mkdir(exist_ok=True)
        with open(outdir / f"{self.kegg.organism}.json", "w") as fout:
            json.dump(self.pathways, fout)

    def _load_pathways(self, progress=True, preload_directory=None):
        # This is just loading all pathways once for all
        self.pathways = {}

        if preload_directory:
            # preload is a directory with all pathways in it
            indir = Path(preload_directory)
            logger.info(
                f"Loading pathways from local files in {preload_directory}. Expecting a file named {self.kegg.organism}.json"
            )
            with open(indir / f"{self.kegg.organism}.json", "r") as fin:
                data = json.load(fin)
                self.pathways = data
            logger.info(f"Loaded {len(self.pathways)} pathways.")
        else:  # pragma: no cover  #not tested due to slow call
            for ID in tqdm(self.kegg.pathwayIds, desc="Downloading KEGG pathways"):
                self.pathways[ID.replace("path:", "")] = self.kegg.parse(self.kegg.get(ID))

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
        go = [x["GO"] if isinstance(x, dict) and "GO" in x.keys() else None for x in self.df_pathways.DBLINKS]
        self.df_pathways["GO"] = go
        del self.df_pathways["DBLINKS"]

    def plot_genesets_hist(self, bins=20):
        N = len(self.gene_sets.keys())
        pylab.clf()
        pylab.hist([len(v) for k, v in self.gene_sets.items()], bins=bins, lw=1, ec="k")
        pylab.title("{} gene sets".format(N))
        pylab.xlabel("Gene set sizes")
        pylab.grid(True, alpha=0.5)
        a, b = pylab.xlim()
        pylab.xlim([0, b])

    def compute_enrichment(self, background=None):
        if background is None:
            background = self.background

        self.summary.data["missing_genes"] = {}
        self.summary.data["input_gene_list"] = {}
        self.enrichment = {
            category: self._enrichr(category, background=background) for category in self.gene_lists.keys()
        }
        self.dfs = {
            category: self._get_final_df(self.enrichment[category].results) for category in self.gene_lists.keys()
        }

    def _enrichr(self, category, background, verbose=True):
        if isinstance(category, list):
            gene_list = category
        else:
            self._check_category(category)
            gene_list = self.gene_lists[category]

        logger.info(f"Computing enrichment for category '{category}'. Input gene list of {len(gene_list)} ids")
        self.summary.data["input_gene_list"][category] = len(gene_list)

        if self.mapper is not None:
            missing = [x for x in gene_list if x not in self.mapper.index]
            logger.info(f"Missing genes from mapper dataframe: {len(missing)}")
            self.summary.data["missing_genes"][category] = ",".join(missing)
            gene_list = [x for x in gene_list if x in self.mapper.index]
            identifiers = self.mapper.loc[gene_list]["name"].drop_duplicates().values

            if self.convert_input_gene_to_upper_case:
                identifiers = [x.upper() for x in identifiers if isinstance(x, str)]
            logger.info(f"Mapped gene list of {len(identifiers)} ids")
            gene_list = list(identifiers)

        gs = GSEA(self.gene_sets)

        enr = gs.compute_enrichment(
            gene_list=gene_list,
            verbose=verbose,
            background=background,
        )

        return enr

    def _get_final_df(self, df):
        # takes the df and populate the name and size of the found pathways
        # we also sort by adjusted p-value
        if len(df) == 0:
            return df

        df = df.copy()
        df["name"] = [self.pathways[x]["NAME"] for x in df.Term]
        df = df.sort_values("Adjusted P-value")

        df["Adjusted P-value (-log10)"] = -np.log10(df["Adjusted P-value"])
        df.reset_index(drop=True, inplace=True)

        df = df.sort_values("Adjusted P-value")
        df = df.rename({"Term": "pathway_id"}, axis=1)
        df = df[df.columns]

        df["significative"] = df["Adjusted P-value"] < self.padj_cutoff
        return df

    def barplot(self, category, nmax=15):
        self._check_category(category)
        df = self.dfs[category].query("significative == True").head(nmax)
        df = df.sort_values("Adjusted P-value (-log10)")

        fig = px.bar(df, x="Adjusted P-value (-log10)", y="name", orientation="h")

        # To have all labels displayed
        if len(df) > 20:
            fig.update_layout(height=800)

        return fig

    def barplot_up_and_down(self, nmax=15):
        import plotly.graph_objects as go

        dfup = self.dfs["up"].query("significative == True").head(nmax)
        dfdown = self.dfs["down"].query("significative == True").head(nmax)
        names = sorted(list(set(list(dfdown["name"].values) + list(dfup["name"].values))))
        names = names[::-1]

        U = [dfup.query("name == @name")["size"].values[0] if name in dfup["name"].values else 0 for name in names]
        D = [-dfdown.query("name == @name")["size"].values[0] if name in dfdown["name"].values else 0 for name in names]

        fig = go.Figure()
        fig.add_trace(go.Bar(y=names, x=U, marker_color="blue", orientation="h", name="Up"))
        fig.add_trace(go.Bar(y=names, x=D, base=0, orientation="h", marker_color="red", name="Down"))

        fig.update_layout(barmode="stack", xaxis_title="Number of Genes")

        return fig

    def scatterplot(self, category, nmax=15):
        self._check_category(category)
        df = self.dfs[category].query("significative == True").head(nmax)

        df.sort_values("Odds Ratio", inplace=True)
        df = df.round({"Odds Ratio": 2, "Combined Score": 2})
        fig = px.scatter(
            df,
            x="Odds Ratio",
            y="name",
            color="Combined Score",
            size="size",
            hover_data=["Combined Score", "Adjusted P-value", "Overlap"],
            color_continuous_scale="Viridis",
        )

        # To have all labels displayed
        if len(df) > 20:
            fig.update_layout(height=800)

        return fig

    def _get_summary_pathway(self, pathway_ID, df, l2fc=0, padj=0.05):
        genes = self.df_pathways.loc[pathway_ID]["GENE"]

        df_down = df.query("padj<=@padj and log2FoldChange<=-@l2fc").copy()
        df_up = df.query("padj<=@padj and log2FoldChange>=@l2fc").copy()

        # for the color, we need to find the match between names and the annotation
        # since names is suppose to have been populated with the annotaion, it should be
        # found. no need for sanity checks.
        df_down["Name"] = df_down[self.color_node_with_annotation]
        df_up["Name"] = df_up[self.color_node_with_annotation]

        # in principle, the KEGG pathays field 'GENE' is stored as
        # "GENE":{key: value} and the values are the NAME, a semi column, and a description
        # hence the split on ; herebelow:
        # unfortunately, there are special cases to handle such as vibrio cholera (vc)
        mapper = {}
        if self.kegg.organism.startswith("vc"):  # pragma: no cover
            for k, v in genes.items():
                # the value is just the description. Let us assume that the name is also the ID
                mapper[k] = k
        else:
            for k, v in genes.items():
                mapper[v.split(";")[0]] = k

        # may not be used ?
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
                    identifier = self.mapper.query("name == @name")["name"].drop_duplicates().index[0]
                    identifiers.append(identifier)
                    new_mapper[identifier] = kegg_id
                except:
                    logger.warning("Skipped {}(kegg ID {}). could not find mapping".format(name, kegg_id))
            mapper = new_mapper

        df_down["Name"] = df_down["Name"].str.lower()
        df_up["Name"] = df_up["Name"].str.lower()

        for name, kegg_id in mapper.items():

            summary_names.append(name)
            summary_keggids.append(kegg_id)

            if name.lower() in {x.lower() for x in df_down.Name.fillna("undefined")}:
                padj = -pylab.log10(df_down.query("Name==@name.lower()").padj.values[0])
                fc = df_down.query("Name==@name.lower()").log2FoldChange.values[0]
                summary_fcs.append(fc)
                summary_pvalues.append(padj)
                summary_types.append("-")
            elif name.lower() in {x.lower() for x in df_up.Name.fillna("undefined")}:
                padj = -pylab.log10(df_up.query("Name==@name.lower()").padj.values[0])
                summary_pvalues.append(padj)
                fc = df_up.query("Name==@name.lower()").log2FoldChange.values[0]
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
        summary["description"] = [self.pathways[pathway_ID]["GENE"][x] for x in summary.keggid]
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
            # to get a value in the range 0,1 expected by the colormap
            l2fc = (l2fc + 4) / 8

            color = colormap.rgb2hex(*cmap(l2fc)[0:3], normalised=True)

            # this is to reduce URL length until a POST solution is found.
            # all entries with low l2fc below |1| are ignored
            # 0.625 <=> l2fc = 1
            # 0.375 <=> l2fc = -1
            if l2fc >= 0.625 or l2fc <= 0.375:
                colors[kegg_id] = f"{color}"

            # 0.125 <==> l2fc = -3 here color is very dark, text should be in another color
            if l2fc < 0.125 or l2fc > 0.875:
                colors[kegg_id] = f"{color},white"
            # else:  #FIXME for now, genes not in the annotation/rnadiff are also ignored
            # to keep URL short
            #    colors[kegg_id] = f"{color}"
        return colors

    def save_pathway(self, pathway_ID, df, scale=None, show=False, filename=None):
        summary = self._get_summary_pathway(pathway_ID, df)
        colors = self._get_colors(summary)

        count_up = len(summary.query("type == '+'"))
        count_down = len(summary.query("type == '-'"))
        logger.info(
            f"pathway {pathway_ID} total genes: {len(summary)}. Found {count_down} down-regulated and {count_up} up-regulated"
        )

        url = "https://www.kegg.jp/kegg-bin/show_pathway"

        params = {
            "map": pathway_ID,
            "multi_query": "\r\n".join([f"{k} {v}" for k, v in colors.items()]),
        }

        if filename is None:
            filename = f"{pathway_ID}.png"

        try:
            # let us try the selenium solution. If it works, we overwrite the previous image
            urlsel = url + f"?map={pathway_ID}&multi_query="
            # print(urlsel)
            # print("\\r\\n".join([f"{k} {v}" for k, v in colors.items()]))
            urlsel += urllib.parse.quote("\r\n".join([f"{k} {v}" for k, v in colors.items()]))
            if len(urlsel) >= 2048:
                logger.warning(f"EXCESS. consider increasing l2fc: {len(urlsel)}, {len(colors)}")
            self._download_kegg_image(urlsel, pathway_ID)
            # in theory one should have an image called PATHWAYID_date_ID.png
            # we need to renae it
            filenames = glob.glob(f"{pathway_ID}@2x_*png")
            if len(filenames) > 1:
                logger.warning(
                    f"expected one PNG file. found several KEGG PNG files: {filenames}. please remove them. Moving the first one found"
                )
            if len(filenames) == 0:
                logger.warning("no files were saved using headless call. roll back to requests ")
            # filename is defined above by the user or as pathway_ID.png
            shutil.move(filenames[0], filename)
        except Exception as err:  # pragma: no cover
            print(err)

            # color do not work anymore since Aug 2023 but this allows us to get the image
            # in case the new implementation (above) fails
            self.params = params
            html_page = requests.post(url, data=params)
            self.tmp = html_page
            html_page = html_page.content.decode()

            try:
                links_to_png = [x for x in html_page.split() if "png" in x and x.startswith("src")]
                link_to_png = links_to_png[0].replace("src=", "").replace('"', "")
                r = requests.get(f"https://www.kegg.jp/{link_to_png}")
                with open(filename, "wb") as fout:
                    fout.write(r.content)
            except IndexError:
                logger.warning(f"could not create image for {pathway_ID}")

        return summary

    def _download_kegg_image(self, url, pathwayID):
        from selenium import webdriver
        from selenium.webdriver.chrome.service import Service as ChromeService
        from selenium.webdriver.common.action_chains import ActionChains
        from selenium.webdriver.common.by import By
        from selenium.webdriver.support import expected_conditions as EC
        from selenium.webdriver.support.ui import WebDriverWait

        # Define the path to the Chrome WebDriver executable
        # Drivers can be found here : https://pypi.org/project/selenium/
        home = expanduser("~")
        webdriver_path = f"{home}/.config/sequana/chromedriver"
        if os.path.exists(webdriver_path) is False:
            logger.critical(
                "To annotate KEGG image, you currently need chrome and a valid/compatible selenium driver for chroe. Please see https://pypi.org/project/selenium"
            )

        # Set up the WebDriver with options
        chrome_options = webdriver.ChromeOptions()
        chrome_options.binary_location = "/usr/bin/google-chrome"  # Path to your Chrome executable
        chrome_options.add_argument("--headless")  # Optional: run Chrome in headless mode (no GUI)

        # Create a ChromeService object with the WebDriver path
        chrome_service = ChromeService(executable_path=webdriver_path)

        logger.info(f"Creating driver to access KEGG website for {pathwayID}")
        # Create a WebDriver instance
        driver = webdriver.Chrome(service=chrome_service, options=chrome_options)

        # Navigate to the web page. Wait 10s for the page to be updated with the requests
        # somehow the WebDriverWait is not enough. Requests larger than 2048 characters will fail
        # here, FIXME with a javascript of post requests but not easy.
        driver.get(url)
        time.sleep(10)

        # Find and click the menu tab 'download-menu'
        download_menu = WebDriverWait(driver, 10).until(EC.presence_of_element_located((By.ID, "downloadMenu")))

        # Use ActionChains to hover over the "Download" menu item
        actions = ActionChains(driver)
        actions.move_to_element(download_menu).perform()

        logger.info("Downloading image\n")
        image_link = WebDriverWait(driver, 10).until(EC.element_to_be_clickable((By.LINK_TEXT, "Image (png) file 2x")))
        actions.move_to_element(image_link).click().perform()

        # Wait for the download to complete (you can use a more robust approach to wait for the file to be ready)
        time.sleep(5)

        # Close the WebDriver
        driver.quit()

    def save_significant_pathways(self, category, nmax=20, tag="", outdir="."):  # pragma: no cover
        """category should be up, down or all"""

        # select the relevant pathways
        df = self.dfs[category].query("significative == True")
        logger.warning("Found {} pathways to save".format(len(df)))
        if len(df) > nmax:
            logger.warning("Restricted pathways to {}".format(nmax))
            df = df.head(nmax)

        logger.info("saving {} deregulated pathways".format(len(df)))

        summaries = {}

        for ID in df["pathway_id"]:
            summary = self.save_pathway(ID, df, filename=(Path(outdir) / f"{ID}_{category}.png"))
            summaries[ID] = summary

        return summaries

    def find_pathways_by_gene(self, gene_name):
        """Returns pathways that contain the gene name

        ::

            ke.find_pathways_by_gene("ysgA")

        """

        # we are going to search for the gene name and lower case onyl
        candidates = [gene_name, gene_name.lower()]

        found_paths = []

        for key in self.pathways.keys():
            if "GENE" in self.pathways[key]:
                genes = set([x.split(";")[0] for x in self.pathways[key]["GENE"].values()])
                for candidate in candidates:
                    if candidate in genes:
                        found_paths.append(key)
        return list(set(found_paths))

    def save_project(self, tag, outdir="."):
        """Save tables and visualisations of the complete enrichmment analysis."""
        outdir = Path(outdir)
        outdir.mkdir(parents=True, exist_ok=True)

        for category in self.gene_lists.keys():
            common_out = Path(f"{tag}_kegg_gsea_{category}_degs")

            df = self.dfs[category]

            if not df.empty:
                df.to_csv(outdir / (common_out.name + ".csv"))
                df.query("significative == True").to_csv(outdir / (common_out.name + "_significant.csv"))

                self.barplot(category).write_image(outdir / (common_out.name + "_barplot.png"))
                self.scatterplot(category).write_image(outdir / (common_out.name + "_scatterplot.png"))

            # In case of no enrichment results, create empty files stating so
            else:
                (outdir / (common_out.name + "_NO_RESULTS")).touch()
