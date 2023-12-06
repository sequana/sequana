#  This file is part of Sequana software
#
#  Copyright (c) 2016-2020 - Sequana Development Team
#
#  Distributed under the terms of the 3-clause BSD license.
#  The full license is in the LICENSE file, distributed with this software.
#
#  website: https://github.com/sequana/sequana
#  documentation: http://sequana.readthedocs.io
#
##############################################################################
import os
import sys
from pathlib import Path

import colorlog
import pandas as pd
import rich_click as click

from sequana.modules_report import ModuleKEGGEnrichment
from sequana.rnadiff import RNADiffResults
from sequana.scripts.utils import CONTEXT_SETTINGS, common_logger
from sequana.utils import config

logger = colorlog.getLogger(__name__)


@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument("rnadiff_directory", type=click.Path(exists=True, file_okay=False, dir_okay=True), nargs=1)
@click.option(
    "--annotation-attribute",
    type=click.STRING,
    default="Name",
    show_default=True,
    help="a valid attribute to be used to map on KEGG database",
)
@click.option(
    "--kegg-name",
    type=click.STRING,
    default=None,
    required=True,
    help="a valid KEGG name (hsa for human, mmu for mus musculus); See the taxonomy command to retrieve other names",
)
@click.option(
    "--log2-foldchange-cutoff",
    type=click.FLOAT,
    default=1,
    show_default=True,
    help="remove events with absolute log2 fold change below this value",
)
@click.option(
    "--padj-cutoff",
    type=click.FLOAT,
    default=0.05,
    show_default=True,
    help="remove events with pvalue above this value default (0.05).",
)
@click.option(
    "--biomart",
    type=click.STRING,
    default=None,
    help="""you may need a biomart mapping of your identifier for the kegg
pathways analysis. If you do not have this file, you can use 'sequana biomart'
command""",
)
@click.option(
    "--plot-linearx",
    type=click.BOOL,
    default=False,
    is_flag=True,
    help="""Default is log2 fold enrichment in the plots. use this to use linear scale""",
)
@click.option(
    "--kegg-pathways-directory",
    type=click.Path(),
    default=None,
    help="""a place where to find the pathways for each organism""",
)
@click.option(
    "--comparison",
    type=click.STRING,
    default=None,
    help="""By default analyses all comparisons found in input file. You may set one specifically with this argument.""",
)
@click.option(
    "--max-pathways",
    type=click.INT,
    default=40,
    show_default=True,
    help="""Max number of pathways to show (most enriched)""",
)
@click.option(
    "--kegg-background",
    type=click.INT,
    default=None,
    help="""a background for kegg enrichment. If None, set to the number of genes
used in the differential analysis (input file rnadiff.csv).""",
)
@click.option(
    "--condition",
    type=click.STRING,
    default="condition",
    help="""The name of the column used in the design file to define groups.""",
)
@click.option("--output-directory", default="enrichment_kegg")
@common_logger
def enrichment_kegg(**kwargs):
    """Create a HTML report showing KEGG enriched pathways
    \b

    Example for the enrichment module:

        sequana enrichment-kegg rnadiff_output_dir --log2-foldchange-cutoff 2

    The KEGG pathways are loaded and it may take time. Once done, they are saved
    in kegg_pathways/organism and be loaded next time:
    \b

        sequana enrichment-kegg rnadiff_output_dir --log2-foldchange-cutoff 2 \\
            --kegg-name lbi --annotation-attribute gene_name


    """
    logger.setLevel(kwargs["logger"])

    keggname = kwargs["kegg_name"]
    params = {
        "padj": kwargs["padj_cutoff"],
        "log2_fc": kwargs["log2_foldchange_cutoff"],
        "mapper": kwargs["biomart"],
        "nmax": kwargs["max_pathways"],
        "kegg_background": kwargs["kegg_background"],
        "preload_directory": kwargs["kegg_pathways_directory"],
        "plot_logx": not kwargs["plot_linearx"],
        "color_node_with_annotation": kwargs["annotation_attribute"],
    }
    filename = kwargs["biomart"]
    if filename and os.path.exists(filename) is False:
        logger.error("{} does not exists".format(filename))
        sys.exit(1)
    filename = kwargs["kegg_pathways_directory"]
    if filename and os.path.exists(filename) is False:
        logger.error("{} does not exists".format(filename))
        sys.exit(1)

    logger.info(f"Reading RNAdiff results from {kwargs['rnadiff_directory']}")
    rnadiff = RNADiffResults(kwargs["rnadiff_directory"], condition=kwargs["condition"])

    # Let us extract the gene names that were used in the
    # differential analysis
    annot_col = kwargs["annotation_attribute"]
    used_genes = rnadiff.annotation[annot_col].loc[rnadiff.df.index].values

    logger.info(f"{len(used_genes)} genes were analysed in the differential analysis")

    # now that we have loaded all results from a rnadiff analysis, let us
    # perform the enrichment for each comparison found in the file
    padj = params["padj"]
    log2fc = params["log2_fc"]

    # setting these attributes set the gene list with log2fc and padj filter
    rnadiff._log2_fc = log2fc
    rnadiff._alpha = padj
    gene_lists = rnadiff.get_gene_lists(
        annot_col=annot_col, Nmax=kwargs.get("max_genes", 1000000), dropna=True
    )  # no filter on number of genes

    output_directory = kwargs["output_directory"]

    if kwargs["comparison"]:
        to_remove = [x for x in gene_lists.keys() if x != kwargs["comparison"]]
        for x in to_remove:
            del gene_lists[x]

    for compa, gene_dict in gene_lists.items():
        if keggname.startswith("vc"):
            logger.warning("Vibrio Cholera annotation old_locus_tag processed to be compatible with KEGG/Sequana")
            if kwargs["annotation_attribute"] == "old_locus_tag":

                def get_kegg_id(x):
                    try:
                        a, b = x.split(",")
                        if a[0:3] == "VC_":
                            return a
                        else:
                            return b
                    except:
                        return x

                for cat in ["up", "down", "all"]:
                    gene_dict[cat] = [get_kegg_id(x) for x in gene_dict[cat]]
                old_locus_tag = rnadiff.annotation["old_locus_tag"].values
                rnadiff.annotation["old_locus_tag"] = [get_kegg_id(x) for x in old_locus_tag]

        config.output_dir = f"{output_directory}/{compa}"
        os.makedirs(f"{output_directory}", exist_ok=True)

        # we define the data and its annotation that will be used by the KEGG
        # enrichment. No need to apply any filter, we pass the entire data set
        # so that even small fold change can be shown
        df = rnadiff.comparisons[compa].df.copy()
        df = pd.concat([df, rnadiff.annotation.loc[df.index].copy()], axis=1)
        df.reset_index(inplace=True)

        m = ModuleKEGGEnrichment(
            gene_dict,
            keggname,
            df,
            enrichment_params=params,
            command=" ".join(["sequana"] + sys.argv[1:]),
            used_genes=used_genes,
        )

    p = Path(f"{output_directory}/.sequana")
    p.mkdir(exist_ok=True)
    with open(p / "info.txt", "w") as fout:
        command = " ".join(["sequana"] + sys.argv[1:])
        fout.write(command)

    m.ke.save_pathways(f"{output_directory}/.sequana")
