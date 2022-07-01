#  This file is part of Sequana software
#
#  Copyright (c) 2016-2022 - Sequana Development Team
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
import json
import itertools

import click
import colorlog

from sequana.enrichment.uniprot_enrichment import UniprotEnrichment
from sequana.modules_report.uniprot_enrichment import ModuleUniprotEnrichment
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
    "--taxon",
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
@click.option("--output-directory", show_default=True, default="enrichment_uniprot")
@common_logger
def enrichment_uniprot(**kwargs):
    """GO enrichment using Uniprot and HTML report creation

    \b
    The input should be the output of the RNADiff pipeline

    Example::

        sequana enrichment-uniprot rnadiff.csv --taxon 10090
            --log2-foldchange-cutoff 2 

        sequana enrichment rnadiff/rnadiff.csv
            --taxon 189518 \
            --log2-foldchange-cutoff 2 
            --ontologies MF 

    \b
    Valid ontologies are: MF, BP, CC


    """
    valid = {"MF", "BP",  "CC"}

    ontologies = eval(kwargs["ontologies"])
    for ontology in ontologies:
        if ontology not in valid:
            logger.error(f"Provided incorrect ontology ({ontology}). Must be in {valid}")
            sys.exit(1)

    logger.setLevel(kwargs["logger"])

    taxon = kwargs["taxon"]
    if taxon == 0:
        logger.error("You must provide a taxon with --taxon")
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
    os.makedirs(output_directory, exist_ok=True)

    for compa, gene_dict in gene_lists.items():
        config.output_dir = f"{output_directory}/{compa}"
        os.makedirs(f"{output_directory}/{compa}", exist_ok=True)

        ue = UniprotEnrichment(gene_dict, 235443)
        ue.compute_enrichment()

        for x in itertools.product(['up', 'down', 'all'], ['CC', 'BP', 'MF']):
            mode, ont = x
            df = ue.plot_go_terms(mode, ontologies=[ont])
            if df is not None and len(df):
                from pylab import savefig
                savefig(f"{output_directory}/{compa}/plot_{mode}_{ont}.png", dpi=200)
                ue.save_chart(df, f'{output_directory}/{compa}/chart_{mode}_{ont}.png')
                df.to_csv(f'{output_directory}/{compa}/{mode}_{ont}.csv', sep=",")


        stats = {"taxon": taxon}
        stats["taxon_name"] = ue.taxon_name
        stats["up"] = len(gene_dict['up'])
        stats["down"] = len(gene_dict['down'])
        stats["all"] = len(gene_dict['all'])
        stats["name"] = compa
        with open(f"{output_directory}/{compa}/summary.json", "w") as fout:
            json.dump(stats, fout, indent=True)


        ModuleUniprotEnrichment(
            gene_dict,
            stats,
            enrichment_params=params,
            command=" ".join(["sequana"] + sys.argv[1:]),
             ontologies=["CC", "BP", "MF"],
        )
