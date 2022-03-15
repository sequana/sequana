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

from sequana.gff3 import GFF3

from .utils import CONTEXT_SETTINGS, common_logger


logger = colorlog.getLogger(__name__)


@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument("input", type=click.Path(exists=True))
@click.argument("output")
@click.option(
    "--features",
    type=click.Path(),
    default="gene",
    help="""list of features to be  extracted""",
)
@common_logger
def gff_to_light_gff(**kwargs):
    """Extract the feature of interest in the input GFF to create a light version

    sequana gff-to-light-gff input.gff output.gff --features gene,exon

    """
    logger.setLevel(kwargs["logger"])
    filename = kwargs["input"]
    assert filename.endswith(".gff") or filename.endswith(".gff3")

    g = GFF3(filename)
    g.read_and_save_selected_features(kwargs["output"], features=kwargs["features"])
