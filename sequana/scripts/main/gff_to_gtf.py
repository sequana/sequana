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

import click
import colorlog

from sequana.gff3 import GFF3

from .utils import CONTEXT_SETTINGS, common_logger


logger = colorlog.getLogger(__name__)


@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument("gff_filename", type=click.Path(exists=True))
@common_logger
def gff_to_gtf(**kwargs):
    """Convert a GFF file into GTF

    This is experimental convertion. Use with care.

    """
    filename = kwargs["gff_filename"]
    assert filename.endswith(".gff") or filename.endswith(".gff3")

    g = GFF3(filename)
    if filename.endswith(".gff"):
        g.to_gtf(os.path.basename(filename).replace(".gff", ".gtf"))
    elif filename.endswith(".gff3"):
        g.to_gtf(os.path.basename(filename).replace(".gff3", ".gtf"))
