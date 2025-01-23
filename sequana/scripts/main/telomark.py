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


"""click.rich_click.OPTION_GROUPS = {
    "telomark": [
        {
            "name": "mosdepth",
            "options": [
                "--window-size",
                "--fast",
            ],
        },
    ],
}
"""


@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument("fasta-file", type=click.STRING)
@click.option(
    "--chromosomes",
    type=click.STRING,
    help="""list of chromosomes to restrict to """,
)
@click.option(
    "--tag",
    type=click.STRING,
    default="telomark",
    help=""" """,
)
@click.option(
    "--chunk-size",
    type=click.INT,
    default=50000,
    help=""" chunk at beginning and end to look at """,
)
@click.option(
    "--peak-height",
    type=click.INT,
    default=20,
    help=""" chunk at beginning and end to look at """,
)
@click.option(
    "--peak-width",
    type=click.INT,
    default=50,
    help=""" chunk at beginning and end to look at """,
)
# @click.option("--exclude-chromosomes", type=click.STRING, default="", help="""list of chromosomes to exclude separated by commas""")
@common_logger
def telomark(**kwargs):
    """ """
    from sequana.telomere import Telomere

    telo = Telomere(kwargs["fasta_file"], peak_height=kwargs["peak_height"], peak_width=kwargs["peak_width"])

    tag = kwargs["tag"]

    if kwargs["chromosomes"]:
        chromosomes = kwargs["chromosomes"].split(",")
        results = telo.run(names=chromosomes, tag=tag, Nmax=kwargs["chunk_size"] * 2)
    else:
        results = telo.run(tag=tag, Nmax=kwargs["chunk_size"] * 2)

    results.to_csv(f"sequana.telomark.{tag}.csv", index=False)

    import pandas as pd
    from pylab import clf, colorbar, imshow, savefig, xlabel, xticks, ylabel, yticks

    clf()
    df = pd.read_csv(f"sequana.telomark.{tag}.csv")
    imshow(
        df[["5to3_LHS_position", "5to3_RHS_position", "3to5_LHS_position", "3to5_RHS_position"]] > 0,
        aspect="auto",
        cmap="YlGn",
    )
    ylabel("chromosome", fontsize=16)

    yticks(range(0, len(df)), df["name"].values, fontsize=8)
    xlabel("telomere", fontsize=16)
    xticks([0, 1, 2, 3], ["LHS", "reversed RHS", "reversed LHS", "RHS"])
    savefig(f"sequana.telomark.summary_binary.{tag}.png", dpi=200)

    clf()
    imshow(
        df[["5to3_LHS_length", "5to3_RHS_length", "3to5_LHS_length", "3to5_RHS_length"]],
        aspect="auto",
        vmin=0,
        vmax=5000,
        cmap="YlGn",
    )
    colorbar()
    yticks(range(0, len(df)), df["name"].values, fontsize=8)
    xlabel("telomere", fontsize=16)
    ylabel("chromosome", fontsize=16)
    xticks([0, 1, 2, 3], ["LHS", "reversed RHS", "reversed LHS", "RHS"])
    savefig(f"sequana.telomark.summary.{tag}.png", dpi=200)
