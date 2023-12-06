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

from sequana.scripts.utils import CONTEXT_SETTINGS

logger = colorlog.getLogger(__name__)


@click.command(context_settings=CONTEXT_SETTINGS)
@click.option("-i", "--input", "infile", help="input FASTA file")
@click.option(
    "-o",
    "--output",
    "outdir",
    default="G4output",
    help="output directory",
)
@click.option("--window", type=click.INT, default=20, show_default=True)
@click.option("--score", type=click.FLOAT, default=1, show_default=True)
def g4hunter(**kwargs):
    """Based on G4Hunter

    takes into account G-richness and G-skewness of a given
    sequence and gives a quadruplexpropensity score as output.')

    """
    from sequana.G4hunter import G4Hunter, G4HunterReader

    G4 = G4Hunter(kwargs["infile"], window=kwargs["window"], score=kwargs["score"])
    G4.run(kwargs["outdir"])
