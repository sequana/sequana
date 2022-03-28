#
#  This file is part of Sequana software
#
#  Copyright (c) 2016-2022 - Sequana Development Team
#
#  Distributed under the terms of the 3-clause BSD license.
#  The full license is in the LICENSE file, distributed with this software.
#
#  website: https://github.com/sequana/sequana
#  documentation: http://sequana.readthedocs.io
#
##############################################################################
import os

from snakemake import shell as shellcmd

import click
import colorlog

from sequana import FastQ
from sequana import FastA

from .utils import CONTEXT_SETTINGS


logger = colorlog.getLogger(__name__)


@click.command(context_settings=CONTEXT_SETTINGS)
@click.option("-1", "--file1", show_default=True, required=True, help="R1 Fastq file ; zipped file expected")
@click.option("-2", "--file2", show_default=True, required=False, help="R2 Fastq file")
@click.option("-r", "--reference", show_default=True, required=True, help="Reference where to map data")
@click.option("-t", "--threads", show_default=True, default=1, help="number of threads to use")
@click.option("-p", "--pacbio", show_default=True, default=False, help="", is_flag=True)
@click.option("-s", "--use-sambamba", is_flag=True, help="use sambamba instad of samtools for sorting")
def mapping(**kwargs):
    """map FastQ data onto a reference

This is a simple mapper for quick test. The commands are as follows:

# Indexing
bwa index REFERENCE\
samtools faidx REFERENCE

# mapping
bwa mem -t 4 -R @RG\\tID:1\\tSM:1\\tPL:illumina -T 30 REFERENCE FASTQ_FILES  | samtools 
view -Sbh -> REFERENCE.bam

samtools sort -o REFERENCE.sorted.bam  REFERENCE.bam 
    """
    reference = kwargs["reference"]
    file1 = kwargs["file1"]
    file2 = kwargs["file2"]
    threads = kwargs["threads"]

    if file1 and file2:
        fastq = f"{file1} {file2}"
    elif file1 and not file2:
        fastq = f"{file1}"
    elif file1 is None:
        raise ValueError("--file1 must be used")

    S = sum([len(this['sequence']) for this in FastQ(file1)])

    if file2:
        S += sum([len(this['sequence']) for this in FastQ(file2)])

    ref = FastA(reference)
    coverage = float(S) / len(ref.sequences[0])
    print("Theoretical Depth of Coverage : %s" % coverage)

    # indexing
    shellcmd(f"bwa index {reference} ")
    cmd = f"samtools faidx {reference} "

    # mapping
    cmd = "bwa mem -M "  # mark shorter split read as secondary; -M is not compulsary but recommended
    if kwargs["pacbio"]:
        cmd += "-x pacbio "
    cmd += f" -t {threads} " + r"-R @RG\\tID:1\\tSM:1\\tPL:illumina -T 30 " + f"{reference} {fastq}  "

    # Samtools options:
    #   S:ignore input format
    #   h:include header
    #   b:bam output
    if kwargs["use_sambamba"] is False:
        cmd += "| samtools view -Sbh | "
        # sorting BAM
        cmd += f"samtools sort -@ {threads} -o {reference}.sorted.bam -"
        shellcmd(cmd)
    else:
        # FIXME use sambamba for the view as well
        cmd += f"| samtools view -Sbu - | sambamba sort /dev/stdin -o {reference}.sorted.bam -t {threads}  --tmpdir=./tmp  "
        shellcmd(cmd)
