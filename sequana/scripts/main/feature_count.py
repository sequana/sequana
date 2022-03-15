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

from .utils import CONTEXT_SETTINGS, common_logger


logger = colorlog.getLogger(__name__)


@click.command(context_settings=CONTEXT_SETTINGS)
@click.option(
    "--pattern",
    help="The pattern of the feature counts files to merge",
    show_default=True,
    default="*feature.out",
)
@click.option(
    "--output",
    help="The output filename where to save the merged counts",
    show_default=True,
    default="all_features.out",
)
@common_logger
def feature_counts(**kwargs):
    """Merge several feature counts files into one file"""
    from sequana.featurecounts import FeatureCountMerger

    fcm = FeatureCountMerger(kwargs["pattern"])
    fcm.to_tsv(output_filename=kwargs["output"])
