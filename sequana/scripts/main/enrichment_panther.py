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

import click
import colorlog

from sequana.modules_report import ModulePantherEnrichment
from sequana.rnadiff import RNADiffResults
from sequana.utils import config

from .utils import CONTEXT_SETTINGS, common_logger, OptionEatAll


logger = colorlog.getLogger(__name__)


# This will be a complex command to provide HTML summary page for
# input files (e.g. bam), or results from pipelines. For each module,
# we should have corresponding option that starts with the module's name
# This can also takes as input various types of data (e.g. FastA)
@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument("name", type=click.Path(exists=True), nargs=1)
@click.option(
    "--annotation-attribute",
    type=click.STRING,
    # required=True,
    default="index",
    show_default=True,
    help="a valid taxon identifiers",
)
@click.option(
    "--panther-taxon",
    type=click.INT,
    required=True,
    help="a valid taxon identifiers",
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
    help="remove events with pvalue abobe this value default (0.05).",
)
@click.option(
    "--plot-linearx",
    type=click.BOOL,
    default=False,
    is_flag=True,
    show_default=True,
    help="""Default is log2 fold enrichment in the plots. use this to use linear scale""",
)
@click.option(
    "--compute-levels/--no-compute-levels",
    default=True,
    help="""Compute the levels of each go term, set --no-compute-levels to skip this step""",
)
@click.option(
    "--max-genes",
    type=click.INT,
    default=2500,
    show_default=True,
    help="""Maximum number of genes (up or down) to use in PantherDB.""",
)
@click.option(
    "--ontologies",
    default=("MF", "BP", "CC"),
    help="""Provide the ontologies to be included in the analysis and HTML report.
Valid choices are: from MF, BP, CC, SLIM_MF, SLIM_BP, SLIM_CC, PROTEIN,
PANTHER_PATHWAY, REACTOME_PATHWAY""",
    cls=OptionEatAll,
    show_default=True,
)
@click.option(
    "--max-enriched-go-terms",
    type=click.INT,
    default=40,
    show_default=True,
    help="""Max number of enriched go terms to show in the plots (most
enriched). All enriched GO terms are stored in tables""",
)
@click.option("--output-directory", show_default=True, default="enrichment_panther")
@common_logger
def enrichment_panther(**kwargs):
    """Create a HTML report for various sequana out

    \b
    * enrichment: the output of RNADiff pipeline

    Example for the enrichment module:

        sequana enrichment-panther rnadiff.csv --panther-taxon 10090
            --log2-foldchange-cutoff 2 

        sequana enrichment rnadiff/rnadiff.csv
            --panther-taxon 189518 \
            --log2-foldchange-cutoff 2 
            --ontologies MF SLIM_MF

    \b
    Valid ontologies are: MF, BP, CC, SLIM_MF, SLIM_BP, SLIM_CC, 
    PROTEIN, "PANTHER_PATHWAY", "REACTOME_PATHWAY"


    """
    valid = [
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

    ontologies = eval(kwargs["ontologies"])
    for ontology in ontologies:
        if ontology not in valid:
            logger.error(f"Provided incorrect ontology ({ontology}). Must be in {valid}")
            sys.exit(1)

    logger.setLevel(kwargs["logger"])

    taxon = kwargs["panther_taxon"]
    if taxon == 0:
        logger.error("You must provide a taxon with --panther-taxon")
        sys.exit(1)

    params = {
        "padj": kwargs["padj_cutoff"],
        "log2_fc": kwargs["log2_foldchange_cutoff"],
        "max_entries": kwargs["max_genes"],
        "nmax": kwargs["max_enriched_go_terms"],
        "plot_logx": not kwargs["plot_linearx"],
        "plot_compute_levels": kwargs["compute_levels"],
    }

    logger.info(f"Reading RNAdiff results from {kwargs['name']}")
    dirpath = os.path.dirname(os.path.abspath(kwargs["name"]))
    rnadiff = RNADiffResults(dirpath, index_col=0, header=[0, 1])

    # now that we have loaded all results from a rnadiff analysis, let us
    # perform the enrichment for each comparison found in the file
    annot_col = kwargs.get("annotation_attribute", "index")

    logger.info(f"Using the annotation column '{annot_col}'")

    # setting these attributes set the gene list with log2fc and padj filter
    rnadiff._log2_fc = params["log2_fc"]
    rnadiff._alpha = params["padj"]
    gene_lists = rnadiff.get_gene_lists(annot_col=annot_col, Nmax=kwargs.get("max_genes", None), dropna=True)

    output_directory = kwargs["output_directory"]
    for compa, gene_dict in gene_lists.items():
        config.output_dir = f"{output_directory}/{compa}"
        os.makedirs(f"{output_directory}", exist_ok=True)

        # for now, let us keep the 'all' category
        # del gene_dict["all"]

        ModulePantherEnrichment(
            gene_dict,
            taxon,
            enrichment_params=params,
            command=" ".join(["sequana"] + sys.argv[1:]),
            ontologies=ontologies,
        )
