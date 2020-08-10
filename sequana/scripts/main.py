# -*- coding: utf-8 -*-

import os
import glob
import json
import sys
import colorlog
import textwrap
import subprocess
import click
from sequana import version
import pathlib

__all__ = ["main"]

from sequana import logger
logger.level = "INFO"


CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])



@click.group(context_settings=CONTEXT_SETTINGS)
@click.version_option(version=version)
def main():
    """This is the main help"""
    pass


@main.command()
@click.argument('filename', type=click.STRING)
@click.option("--count-reads", is_flag=True)
@click.option("-o", "--output", 
    help="filename where to save results. to be used with --head, --tail")
@click.option("--head", type=click.INT, 
    help='number of reads to extract from the head')
@click.option("--tail", type=click.INT, 
    help="number of reads to extract from the tail")
def fastq(**kwargs):
    """Set of useful utilities for FastQ manipulation.

    Input file can be gzipped or not. The --output-file

    """
    from sequana.fastq import FastQ

    # could be simplified calling count_reads only once
    if kwargs['count_reads']:
        import glob
        for filename in glob.glob(kwargs["filename"]):
            f = FastQ(filename)
            Nreads = f.count_reads()
            Nlines = Nreads * 4
            print(f"Number of reads in {filename}: {Nreads}")
            print(f"Number of lines in {filename}: {Nlines}")
    elif kwargs['head']:
        f = FastQ(kwargs['filename'])
        assert kwargs['output']
        N = kwargs['head'] * 4 
        f.extract_head(N=N, output_filename=kwargs['output'])
    elif kwargs['tail']:
        raise NotImplementedError
    else:
        print("Use one of the commands")




@main.command()
@click.argument('name', type=click.STRING)
def develop(**kwargs):
    """

    """
    name = kwargs['name']


if __name__ == "__main__": #pragma: no cover
    main()
