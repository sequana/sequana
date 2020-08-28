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


import pkg_resources
pipelines = [item.key for item in pkg_resources.working_set if item.key.startswith("sequana")]
if len(pipelines): 
    version +="\nThe following pipelines are installed:\n"
for item in pkg_resources.working_set:
    if item.key.startswith("sequana") and item.key != 'sequana':
        version += "\n - {}, version {}".format(item.key, item.version)

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
    elif kwargs['tail']: #pragma: no cover
        raise NotImplementedError
    else:  #pragma: no cover
        print("Use one of the commands")


@main.command()
@click.argument('name', type=click.STRING)
@click.option('--check', is_flag=True)
@click.option('--extract-adapters', is_flag=True)
def samplesheet(**kwargs):
    """Utilities to manipulate sample sheet"""
    name = kwargs['name']
    from sequana.iem import IEM
    iem = IEM(name)
    if kwargs['check']:
        iem.validate()
    elif kwargs["extract_adapters"]:
        iem.to_fasta()

@main.command()
@click.argument("name", type=click.STRING)
@click.option("--module")
@click.option("--enrichment-taxon", type=click.INT)
def summary(**kwargs):
    name = kwargs['name']
    # we will try by monkey patching
    from sequana.rnadiff import RNADiffResults
    from sequana.modules_report.rnadiff import RNAdiffModule
    from sequana.modules_report.bamqc import BAMQCModule

    module = kwargs['module']
    if module:
        if module == "bamqc":
            report = BAMQCModule(name, "bamqc.html")
        elif module == "rnadiff":
            data = RNADiffResults(name)
            report = RNAdiffModule(data)
        elif module == "enrichment":
            from sequana.modules_report.enrichment import Enrichment
            report = Enrichment(data, taxon)
    else: # we try everthing
        found = False
        if found is False:
            try:
                data = RNADiffResults(name)
                report = RNAdiffModule(data)
                found = True
                print("RNADiff module created ")
            except Exception as err: 
                pass

        if found is False:
            try: 
                report = BAMQCModule(name, "bamqc.html")
                found = True
                print("BAMQC module created in bamqc.html")
            except:pass


@main.command()
@click.option("--input", required=True,
    help="The salmon input file.")
@click.option("--output", required=True,
    help="The feature counts output file")
@click.option("--gff", required=True,
    help="A GFF file compatible with your salmon file")
@click.option("--attribute", default="ID",
    help="A valid attribute to be found in the GFF file and salmon input")
def salmon(**kwargs):
    """Convert output of Salmon into a feature counts file """
    from sequana import salmon
    salmon_input = kwargs['input']
    output = kwargs["output"]
    if os.path.exists(salmon_input) is False:
        logger.critical("Input file does not exists ({})".format(salmon_input))
    gff = kwargs["gff"]
    attribute = kwargs['attribute']
    s = salmon.Salmon(salmon_input)
    s.save_feature_counts(output, gff, attribute=attribute)


if __name__ == "__main__": #pragma: no cover
    main()

