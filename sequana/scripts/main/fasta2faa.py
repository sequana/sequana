#  This file is part of Sequana software
#
#  Copyright (c) 2026 - Sequana Development Team
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

from sequana.codecs.fasta2faa import fasta2faa as conv
from sequana.scripts.utils import CONTEXT_SETTINGS, common_logger

logger = colorlog.getLogger(__name__)


@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument("input_fasta", type=click.Path(exists=True))
@click.argument("output_faa", type=click.Path())
@click.option("--width", type=click.INT, default=60, help="Line-wrap width for the output protein sequences.")
@common_logger
def fasta2faa(**kwargs):
    """Convert a nucleotide FASTA file to a protein FASTA (.faa) file.

    Naive frame-0 translation of every sequence (no ORF detection, no
    reverse-complement, no alternative reading frames). Stop codons are
    written as '_'; ambiguous or trailing partial codons as 'X'.
    """
    input_fasta = kwargs["input_fasta"]
    output_faa = kwargs["output_faa"]
    width = kwargs["width"]

    count = conv(input_fasta, output_faa, width=width)
    logger.info(f"Translated {count} sequence(s) from {input_fasta} to {output_faa}")
