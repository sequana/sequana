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

from sequana.gtf_fixer import GTFFixer
from sequana.scripts.utils import CONTEXT_SETTINGS

logger = colorlog.getLogger(__name__)


@click.command(context_settings=CONTEXT_SETTINGS)
@click.option("-i", "--input", required=True)
@click.option("-o", "--output", required=True)
def gtf_fixer(**kwargs):
    """Reads GTF and fix known issues (exon and genes uniqueness)"""
    gtf = GTFFixer(kwargs["input"])
    res = gtf.fix(kwargs["output"])
    print(res)
