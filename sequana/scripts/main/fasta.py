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
import glob
import os
import subprocess
import sys

import click
import colorlog

from sequana import FastA

from .utils import CONTEXT_SETTINGS


logger = colorlog.getLogger(__name__)


@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument("filename", type=click.STRING, nargs=-1)
@click.option(
    "-o",
    "--output",
    help="filename where to save results. to be used with --head, --tail",
)
@click.option("--count-sequences", is_flag=True)
@click.option("--head", type=click.INT, help="number of sequences to extract from the head")
@click.option("--merge", is_flag=True, help="merge all compressed input fastq files into a single file")
@click.option("--tail", type=click.INT, help="number of reads to extract from the tail")
@click.option("--explode", is_flag=True, help="Create a fasta file for each sequence found in the original files")
def fasta(**kwargs):
    """Set of useful utilities for FastA manipulation."""
    filenames = kwargs["filename"]
    # users may provide a wildcards such as "A*gz" or list of files.
    if len(filenames) == 1:
        # if existing files or glob, a glob would give the same answer.
        filenames = glob.glob(filenames[0])
    for filename in filenames:
        os.path.exists(filename)

    # could be simplified calling count_reads only once
    if kwargs["count_sequences"]:
        for filename in filenames:
            f = FastA(filename)
            Nreads = len(f)
            print(f"Number of reads in {filename}: {Nreads}")
    elif kwargs["merge"]:
        # merge all input files (assuming gz extension)
        extensions = [filename.split(".")[-1] for filename in filenames]
        if set(extensions) != set(["gz"]):
            raise ValueError("Your input FastA files must be zipped")
        output_filename = kwargs["output"]
        if output_filename is None:
            logger.error("You must use --output filename.gz")
            sys.exit(1)
        if output_filename.endswith(".gz") is False:
            raise ValueError("your output file must end in .gz")
        p1 = subprocess.Popen(["zcat"] + list(filenames), stdout=subprocess.PIPE)
        fout = open(output_filename, "wb")
        subprocess.run(["pigz"], stdin=p1.stdout, stdout=fout)
    elif kwargs["explode"]:
        for filename in filenames:
            f = FastA(filename)
            f.explode()
