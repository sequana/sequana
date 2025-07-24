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
import sys
from pathlib import Path

import colorlog
import rich_click as click
from tqdm import tqdm

from sequana.scripts.common import teardown
from sequana.scripts.utils import CONTEXT_SETTINGS, common_logger
from sequana.utils import config

logger = colorlog.getLogger(__name__)


click.rich_click.OPTION_GROUPS = {
    "somy score": [
        {
            "name": "mosdepth",
            "options": [
                "--window-size",
                "--fast",
            ],
        },
    ],
}


@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument("filename", type=click.STRING)
@click.option(
    "--window-size",
    type=click.INT,
    default=1000,
    show_default=True,
    help="""The reference to test DGE against. If provided, conditions not
            involving the reference are ignored. Otherwise all combinations are
            tested""",
)
@click.option(
    "--fast",
    is_flag=True,
    help="""fast option""",
)
@click.option(
    "--method",
    type=click.Choice(["em", "median", "mean"]),
    default="median",
    show_default=True,
    help="""Method to estimate the main somy (in principle diploid). Median and mean methods simply compute those statistics from the chunks (windows). The EM is more complex, and compute distribution, estimate mixture model and therefore the mean of the diploid distribution based on those estimates""",
)
@click.option(
    "--estimated-diploy-coverage",
    default=None,
    required=False,
    help="""If not provided, data normalisation is based on --method. If provided, this is the estimated coverage""",
)
@click.option(
    "--chromosomes",
    type=click.STRING,
    help="""list of chromosomes to restrict to """,
)
@click.option(
    "--mapq",
    type=click.INT,
    default=0,
    show_default=True,
    help="""list of chromosomes to restrict to """,
)
@click.option(
    "--telomeric-span",
    type=click.INT,
    default=10,
    show_default=True,
    help="""region to ignore in kb. This suppose that contigs are circularised correctly. If not, set to 0""",
)
@click.option(
    "-k",
    type=click.INT,
    default=4,
    show_default=True,
    help="""Model for gaussin mixture (k=4 is suppose to capture di + tri + tetraploidy)""",
)
@click.option(
    "--minimum-depth", type=click.INT, default=10, help="""coverage with depth below this value are removed."""
)
@click.option(
    "--flag",
    type=click.INT,
    default=3844,
    show_default=True,
    help="""3844 means that it removes unmapped reads, but also secondary and supplementary reads.""",
)
@click.option("--threads", type=click.INT, default=4, help="""use 4 threads .""")
@click.option("--exclude-chromosomes", type=click.STRING, default="", help="""list of chromosomes to exclude""")
@common_logger
def somy_score(**kwargs):
    """**Somy score on polyploid (or not)**"""
    import pandas as pd
    from easydev import cmd_exists

    from sequana import logger
    from sequana.modules_report.rnadiff import RNAdiffModule

    logger.setLevel(kwargs["logger"])

    if not cmd_exists("mosdepth"):
        logger.critical(
            """mosdepth not found. You may install it yourself or use damona using 'damona install mosdepth' """
        )
        sys.exit(1)

    exclude_chromosomes = kwargs.get("exclude_chromosomes", "").split(",")

    from sequana.somy import SomyScore

    ss = SomyScore(kwargs["filename"], window_size=kwargs["window_size"])
    if kwargs["chromosomes"]:
        ss.compute_coverage(
            use_median=True,
            fast=True,
            mapq=kwargs["mapq"],
            flag=kwargs["flag"],
            chromosomes=kwargs["chromosomes"].split(","),
            threads=kwargs["threads"],
            exclude_chromosomes=exclude_chromosomes,
        )
    else:
        ss.compute_coverage(
            use_median=True,
            fast=True,
            mapq=kwargs["mapq"],
            flag=kwargs["flag"],
            threads=kwargs["threads"],
            exclude_chromosomes=exclude_chromosomes,
        )
    # save results before filtering so that we can introspect the data later on
    logger.info(f"length input coverage file: {len(ss.df)}. Saving data into data.csv")
    ss.df.to_csv("data.csv", index=False)

    ss.remove_outliers()
    logger.info(f"Removed outliers. length input coverage file: {len(ss.df)}")

    ss.remove_low_depth(kwargs["minimum_depth"])
    logger.info(f"Removed low depth. length input coverage file: {len(ss.df)}")

    flanks = kwargs.get("telomeric_span")
    ss.remove_flanks(remove_flanking_regions_kb=flanks)
    logger.info(f"Removing flanks ({flanks}kb). length input coverage file: {len(ss.df)}")

    from matplotlib import rcParams

    rcParams["figure.figsize"] = (12, 8)

    if kwargs["estimated_diploy_coverage"]:
        ss.boxplot(
            k=kwargs["k"], method=kwargs["method"], hybrid=True, muhat=float(kwargs["estimated_diploy_coverage"])
        )
    else:
        ss.boxplot(k=kwargs["k"], method=kwargs["method"], hybrid=True)

    print(ss.info)
