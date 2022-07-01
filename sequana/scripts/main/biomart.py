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
import click
import colorlog

from .utils import CONTEXT_SETTINGS, common_logger


logger = colorlog.getLogger(__name__)


@click.command(context_settings=CONTEXT_SETTINGS)
@click.option(
    "--mart",
    default="ENSEMBL_MART_ENSEMBL",
    show_default=True,
    help="A valid mart name",
)
@click.option(
    "--host",
    default=None,
    show_default=True,
    help="A valid mart name such as www.ensembl.org (default) or fungi.ensembl.org",
)
@click.option(
    "--dataset",
    required=True,
    help="A valid dataset name. e.g. mmusculus_gene_ensembl, hsapiens_gene_ensembl",
)
@click.option(
    "--attributes",
    default="ensembl_gene_id,go_id,entrezgene_id,external_gene_name",
    show_default=True,
    help="""A valid set of attributes to look for in the dataset. Multiple
attributes are separeted by a comma (no spaces accepted)""",
)
@click.option(
    "--output",
    default=None,
    help="""by default save results into a CSV file named
    biomart_<dataset>_<YEAR>_<MONTH>_<DAY>.csv""",
)
@common_logger
def biomart(**kwargs):
    """Retrieve information from biomart and save into CSV file

    This command uses BioMart from BioServices to introspect a MART service
    (--mart) and a specific dataset (default to mmusculus_gene_ensembl). Then,
    for all ensembl IDs, it will fetch the requested attributes (--attributes).
    Finally, it saves the CSV file into an output file (--output). This takes
    about 5-10 minutes to retrieve the data depending on the connection.

    Example:

        sequana biomart --mart mmusculus_gene_ensembl --mart ENSEMBL_MART_ENSEMBL \
            --dataset mmusculus_gene_ensembl \
            --attributes ensembl_gene_id,external_gene_name,go_id \
            --output test.csv

    """
    logger.setLevel(kwargs["logger"])

    mart = kwargs["mart"]
    attributes = kwargs["attributes"]
    dataset = kwargs["dataset"]
    host = kwargs["host"]

    from sequana.enrichment.mart import Mart

    conv = Mart(host=host, dataset=dataset, mart=mart)

    df = conv.query(attributes.split(","))
    conv.save(df, filename=kwargs["output"])
