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
import sys

import click
import colorlog
import shutil

from sequana.ribodesigner import RiboDesigner

from .utils import CONTEXT_SETTINGS


logger = colorlog.getLogger(__name__)


# =====================================================================================
# Ribodepletion custom probes designer
# =====================================================================================
@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument("fasta", type=click.Path(exists=True))
@click.argument("gff", type=click.Path(exists=True))
@click.option("--output-directory", default="out_ribodesigner", type=click.Path(exists=False))
@click.option("--seq-type", default="rRNA", help="The annotation type (column 3 in gff) to target for probes.")
@click.option("--max-n-probes", default=384, type=click.INT, help="The maximum number of probes to design.")
@click.option("--threads", default=4, type=click.INT, help="The number of threads to use for cd-hit-est.")
@click.option("--identity-step", default=0.01, type=click.FLOAT, help="The number of threads to use for cd-hit-est.")
@click.option("--output-image", default=None)
@click.option(
    "--force",
    default=False,
    is_flag=True,
    type=click.BOOL,
    help="If output directory exists, use this option to erase previous results",
)
def ribodesigner(**kwargs):
    """A tool to design custom ribodepletion probes.

    This uses a reference genome (FASTA file) and the corresponding annotation
    (GFF file). CD-HIT-EST should be installed and in your $PATH.
    """

    if not shutil.which("cd-hit-est"):
        logger.error("cd-hit-est not found in PATH.")
        sys.exit(1)

    RiboDesigner(**kwargs).run()
    if kwargs["output_image"]:
        from pylab import savefig
        savefig("/".join([kwargs["output_directory"], kwargs["output_image"]]))
