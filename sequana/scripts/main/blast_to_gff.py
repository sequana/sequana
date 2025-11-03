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

import colorlog
import rich_click as click

from sequana.blast import blast_to_gff as blast2gff
from sequana.scripts.utils import CONTEXT_SETTINGS, common_logger

logger = colorlog.getLogger(__name__)


@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument("input_blast", type=click.Path(exists=True))
@click.argument("output_gff", type=click.Path())
@common_logger
def blast_to_gff(**kwargs):
    """Convert a Blast file (outfmt=6)"""
    blast = kwargs["input_blast"]
    gff = kwargs["output_gff"]
    blast2gff(blast, gff)
