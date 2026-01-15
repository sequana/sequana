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
import glob
import os
import subprocess
import sys

import colorlog
import rich_click as click

from sequana.gff3 import GFF3
from sequana.scripts.utils import CONTEXT_SETTINGS

logger = colorlog.getLogger(__name__)


@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument("filename", type=click.STRING, nargs=-1)
@click.option(
    "-o",
    "--output",
    help="filename where to save results.",
)
@click.option("--add-CDS-and-mRNA", is_flag=True, type=click.BOOL, help="add CDS and mRNA")
@click.option(
    "--gene-id",
    type=click.STRING,
    help="Uses this identifier to get gene names to build new DS and mRNA identifiers",
    default="gene_id",
)
def gff(**kwargs):
    """Set of useful utilities for GFF manipulation."""
    filenames = kwargs["filename"]
    # users may provide a wildcards such as "A*gz" or list of files.
    if len(filenames) == 1:
        # if existing files or glob, a glob would give the same answer.
        filenames = glob.glob(filenames[0])
    for filename in filenames:
        os.path.exists(filename)

    # could be simplified calling count_reads only once
    if kwargs["add_cds_and_mrna"]:
        for filename in filenames:
            g = GFF3(filename)
            g.add_CDS_and_mRNA(output=kwargs["output"], gene_ID=kwargs["gene_id"])
