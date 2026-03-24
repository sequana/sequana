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

import json
import xml.etree.ElementTree as ET
from io import BytesIO
from pathlib import Path

import colorlog
import requests
from tqdm import tqdm

from sequana.enrichment.gsea import GSEA
from sequana.lazy import bioservices, colormap
from sequana.lazy import numpy as np
from sequana.lazy import pandas as pd
from sequana.lazy import plotly_express as px
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

        self.kegg = bioservices.KEGG(cache=True)
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

        if filename is None:
            filename = f"{pathway_ID}.png"

        try:
            self._render_pathway_image(pathway_ID, colors, filename)
        except Exception as err:
            logger.warning(f"Could not create image for {pathway_ID}: {err}")

        return summary

    def _render_pathway_image(self, pathway_id, colors, filename):
        """Download KEGG pathway static image and overlay gene colors using KGML coordinates.

        Uses the KEGG REST API:
        - ``https://rest.kegg.jp/get/{pathway_id}/image`` for the background PNG
        - ``https://rest.kegg.jp/get/{pathway_id}/kgml`` for node bounding boxes

        :param pathway_id: KEGG pathway ID (e.g. ``eco00010``)
        :param colors: dict mapping KEGG gene entry IDs to hex color strings
            (optionally ``"#rrggbb,white"`` for contrasting text, only the hex part is used)
        :param filename: output PNG file path
        """
        from matplotlib import colormaps
        from PIL import Image, ImageDraw, ImageFont

        try:
            cmap = colormaps["PuOr_r"]
        except Exception:
            from matplotlib.cm import get_cmap

            cmap = get_cmap("PuOr_r")

        img_resp = requests.get(f"https://rest.kegg.jp/get/{pathway_id}/image")
        img_resp.raise_for_status()
        img = Image.open(BytesIO(img_resp.content)).convert("RGBA")

        kgml_resp = requests.get(f"https://rest.kegg.jp/get/{pathway_id}/kgml")
        kgml_resp.raise_for_status()
        root = ET.fromstring(kgml_resp.text)

        overlay = Image.new("RGBA", img.size, (0, 0, 0, 0))
        draw = ImageDraw.Draw(overlay)

        try:
            font = ImageFont.load_default(size=9)
        except TypeError:
            font = ImageFont.load_default()

        logger.info(f"colors dict: {len(colors)} entries, sample keys: {list(colors.keys())[:5]}")

        gene_entries = root.findall("entry[@type='gene']")
        logger.info(f"KGML gene entries: {len(gene_entries)}")

        _RED_BORDER = (200, 0, 0, 255)  # red outline for map/pathway-link boxes
        _PURPLE_BORDER = (128, 0, 128, 255)  # purple outline for compound/metabolite circles

        def _draw_node(gx, gy, gw, gh, bg, text_col, label, shape="rectangle", outline=(0, 0, 0, 255), outline_width=1):
            box = (gx - gw // 2, gy - gh // 2, gx + gw // 2, gy + gh // 2)
            if shape == "circle":
                draw.ellipse(box, fill=bg, outline=outline, width=outline_width)
            elif shape == "roundrectangle":
                draw.rounded_rectangle(box, radius=4, fill=bg, outline=outline, width=outline_width)
            else:
                draw.rectangle(box, fill=bg, outline=outline, width=outline_width)
            if label:
                lbbox = draw.textbbox((0, 0), label, font=font)
                lw, lh = lbbox[2] - lbbox[0], lbbox[3] - lbbox[1]
                draw.text((gx - lw // 2, gy - lh // 2), label, fill=text_col, font=font)

        # Pass 1a: paint all gene nodes white so DE colors stand out
        for entry in gene_entries:
            gr = entry.find("graphics")
            if gr is None:
                continue
            _draw_node(
                int(gr.get("x", 0)),
                int(gr.get("y", 0)),
                int(gr.get("width", 46)),
                int(gr.get("height", 17)),
                bg=(255, 255, 255, 255),
                text_col=(0, 0, 0, 255),
                label=gr.get("name", "").split(",")[0].strip(),
            )

        # Pass 1b: paint ortholog boxes and compound circles green to restore visual variety
        for etype in ("ortholog", "compound"):
            for entry in root.findall(f"entry[@type='{etype}']"):
                gr = entry.find("graphics")
                if gr is None:
                    continue
                shape = gr.get("type", "rectangle")
                label = ""  # original text already in the KEGG image; avoid overwriting
                _draw_node(
                    int(gr.get("x", 0)),
                    int(gr.get("y", 0)),
                    int(gr.get("width", 8)),
                    int(gr.get("height", 8)),
                    bg=(0, 0, 0, 0),  # transparent fill
                    text_col=(0, 0, 0, 255),
                    label=label,
                    shape=shape,
                    outline=_PURPLE_BORDER,
                    outline_width=2,
                )

        # Pass 1c: map (pathway link) boxes — red border only, no fill so original text shows through
        for entry in root.findall("entry[@type='map']"):
            gr = entry.find("graphics")
            if gr is None:
                continue
            _draw_node(
                int(gr.get("x", 0)),
                int(gr.get("y", 0)),
                int(gr.get("width", 100)),
                int(gr.get("height", 25)),
                bg=(0, 0, 0, 0),  # transparent fill
                text_col=(0, 0, 0, 255),
                label="",
                shape=gr.get("type", "roundrectangle"),
                outline=_RED_BORDER,
                outline_width=3,
            )

        # Pass 2: color DE genes
        colored = 0
        for entry in gene_entries:
            gene_ids = [gid.split(":")[-1] for gid in entry.get("name", "").split()]
            gr = entry.find("graphics")
            if gr is None:
                continue
            for gid in gene_ids:
                if gid in colors:
                    color_str = colors[gid]
                    parts = color_str.split(",")
                    hex_color = parts[0]
                    text_col = (255, 255, 255, 255) if len(parts) > 1 else (0, 0, 0, 255)
                    cr, cg, cb = int(hex_color[1:3], 16), int(hex_color[3:5], 16), int(hex_color[5:7], 16)
                    _draw_node(
                        int(gr.get("x", 0)),
                        int(gr.get("y", 0)),
                        int(gr.get("width", 46)),
                        int(gr.get("height", 17)),
                        bg=(cr, cg, cb, 255),
                        text_col=text_col,
                        label=gr.get("name", "").split(",")[0].strip(),
                    )
                    colored += 1
                    break

        # Pass 3: draw legend
        self._draw_log2fc_legend(draw, img.size, font, cmap)

        logger.info(f"Colored {colored}/{len(gene_entries)} gene nodes in {pathway_id}")
        Image.alpha_composite(img, overlay).convert("RGB").save(filename)

    def _draw_log2fc_legend(self, draw, img_size, font, cmap):
        """Draw a log2FC color gradient legend in the bottom-right corner of the overlay."""
        from PIL import ImageFont

        # Try to load a TrueType font for crisp rendering; fall back to default bitmap font
        _FONT_PATHS = [
            "/usr/share/fonts/dejavu-sans-fonts/DejaVuSans.ttf",
            "/usr/share/fonts/dejavu/DejaVuSans.ttf",
            "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf",
            "/usr/share/fonts/liberation-sans/LiberationSans-Regular.ttf",
            "/usr/share/fonts/truetype/liberation/LiberationSans-Regular.ttf",
            "/usr/share/fonts/google-noto/NotoSans-Regular.ttf",
        ]
        legend_font = None
        for fp in _FONT_PATHS:
            try:
                legend_font = ImageFont.truetype(fp, 11)
                break
            except (IOError, OSError):
                continue
        if legend_font is None:
            legend_font = font  # bitmap fallback

        bar_w, bar_h = 200, 16
        pad = 10

        sample_bbox = draw.textbbox((0, 0), "0", font=legend_font)
        th = sample_bbox[3] - sample_bbox[1]  # text height

        # Layout: title / gap / bar / gap / tick labels
        content_h = th + 4 + bar_h + 4 + th
        img_w, img_h = img_size
        lx = img_w - bar_w - pad - 10  # left edge of bar
        ly = pad + 8  # top of title — top-right corner

        # Semi-transparent white background box
        draw.rectangle(
            (lx - 8, ly - 5, lx + bar_w + 8, ly + content_h + 5),
            fill=(255, 255, 255, 220),
            outline=(100, 100, 100, 255),
        )

        # Title — use plain "log2FC" to avoid Unicode rendering issues with bitmap fonts
        title = "log2FC"
        tbbox = draw.textbbox((0, 0), title, font=legend_font)
        tw = tbbox[2] - tbbox[0]
        draw.text((lx + bar_w // 2 - tw // 2, ly), title, fill=(40, 40, 40, 255), font=legend_font)

        # Gradient bar
        bar_y = ly + th + 4
        for i in range(bar_w):
            rgba = cmap(i / (bar_w - 1))
            draw.line(
                [(lx + i, bar_y), (lx + i, bar_y + bar_h - 1)],
                fill=(int(rgba[0] * 255), int(rgba[1] * 255), int(rgba[2] * 255), 255),
            )
        draw.rectangle([lx, bar_y, lx + bar_w, bar_y + bar_h], outline=(80, 80, 80, 255))

        # Tick marks and labels below bar
        tick_y = bar_y + bar_h + 4
        for l2fc_val, tick_label in [(-4, "-4"), (-1, "-1"), (0, "0"), (1, "+1"), (4, "+4")]:
            pos = int((l2fc_val + 4) / 8 * (bar_w - 1))
            tx = lx + pos
            draw.line([(tx, bar_y + bar_h), (tx, bar_y + bar_h + 3)], fill=(80, 80, 80, 255))
            lbbox = draw.textbbox((0, 0), tick_label, font=legend_font)
            lw = lbbox[2] - lbbox[0]
            draw.text((tx - lw // 2, tick_y), tick_label, fill=(40, 40, 40, 255), font=legend_font)

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
