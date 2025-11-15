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

from sequana.codecs.embl_to_fasta import embl_to_fasta as conv
from sequana.scripts.utils import CONTEXT_SETTINGS, common_logger

logger = colorlog.getLogger(__name__)


@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument("input_embl", type=click.Path(exists=True))
@click.argument("output_fasta", type=click.Path())
@common_logger
def embl_to_fasta(**kwargs):
    """Convert a emble to fasta"""
    blast = kwargs["input_embl"]
    gff = kwargs["output_fasta"]
    conv(blast, gff)
