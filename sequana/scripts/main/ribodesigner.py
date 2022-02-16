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
from pathlib import Path

import click
import colorlog
from easydev import cmd_exists

from sequana import ribodesigner as rd

from .utils import CONTEXT_SETTINGS


logger = colorlog.getLogger(__name__)


@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument("fasta", type=click.Path(exists=True))
@click.argument("gff", type=click.Path(exists=True))
@click.option("--output-directory", default="out_ribodesigner", type=click.Path(exists=False))
@click.option("--seq-type", default="rRNA", help="The annotation type (column 3 in gff) to target for probes.")
@click.option("--probe-len", default=50, type=click.INT, help="The size of probes in nucleotides.")
@click.option(
    "--inter-probe-space", default=15, type=click.INT, help="The size of the space between probes in nucleotides."
)
@click.option("--seq-id-thres", default=0.80, type=click.FLOAT, help="The similarity threshold for clustering.")
@click.option("--threads", default=4, type=click.INT, help="The number of threads to use for cd-hit-est.")
def ribodesigner(**kwargs):
    """A tool to design custom ribodepletion probes.

    This uses a reference genome (FASTA file) and the corresponding annotation
    (GFF file). CD-HIT-EST should be installed and in your $PATH.
    """
    if not cmd_exists("cd-hit-est"):
        logger.error("cd-hit-est not found in PATH.")
        sys.exit(1)

    outdir = Path(kwargs["output_directory"])
    outdir.mkdir()

    regions_fas = outdir / "regions.fas"
    probes_fas = outdir / "probes.fas"
    clustered_fas = outdir / "clustered_probes.fas"
    clustered_csv = outdir / "clustered_probes.csv"

    rd.get_rna_pos_from_gff(
        kwargs["gff"], outdir / "annot_filtered.gff", kwargs["fasta"], regions_fas, seq_type=kwargs["seq_type"]
    )

    # rd.extract_regions_from_fasta(kwargs["fasta"], df, regions_fas)
    rd.get_probes(regions_fas, probes_fas, probe_len=kwargs["probe_len"], inter_probe_space=kwargs["inter_probe_space"])
    rd.cluster_probes(probes_fas, clustered_fas, seq_id_thres=kwargs["seq_id_thres"], threads=kwargs["threads"])
    rd.fasta_to_csv(clustered_fas, clustered_csv)
