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
from pylab import savefig

from sequana.compare import RNADiffCompare

from .utils import CONTEXT_SETTINGS, common_logger


logger = colorlog.getLogger(__name__)


@click.command(context_settings=CONTEXT_SETTINGS)
@click.option(
    "--file1",
    type=click.Path(),
    default=None,
    required=True,
    help="""The first input RNA-seq table to compare""",
)
@click.option(
    "--file2",
    type=click.Path(),
    default=None,
    required=True,
    help="""The second input RNA-seq table to compare""",
)
@common_logger
def rnaseq_compare(**kwargs):
    """Compare 2 tables created by the 'sequana rnadiff' command"""
    c = RNADiffCompare(kwargs["file1"], kwargs["file2"])
    print(c.r1.summary())
    print(c.r2.summary())
    c.plot_volcano_differences()
    savefig("sequana_rnaseq_compare_volcano.png", dpi=200)
