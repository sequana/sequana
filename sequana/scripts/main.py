# -*- coding: utf-8 -*-

import sys
import os
import glob
import click
from sequana import version
import functools

__all__ = ["main"]

import sequana
from sequana import logger

logger.level = "INFO"

# This can be used by all commands as a simple decorator
def common_logger(func):
    @click.option("--logger", default="INFO",
        type=click.Choice(["INFO", "DEBUG", "WARNING", "CRITICAL", "ERROR"]))
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        return func(*args, **kwargs)
    return wrapper


def get_env_vars(ctx, args, incomplete):
    return [k for k in os.environ.keys() if incomplete in k]

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
def main(**kwargs):
    """\bThis is the main entry point for a set of Sequana applications.

    Pipelines such as sequana_rnaseq, sequana_variant_calling have their own
    application and help.

    In addition, more advanced tools such as sequana_taxonomy or
    sequana_coverage have their own standalone.

    """
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


# This will be a complex command to provide HTML summary page for
# input files (e.g. bam), or results from pipelines. For each module,
# we should have corresponding option that starts with the module's name
@main.command()
@click.argument("name", type=click.Path(exists=True))
@click.option("--module",
    required=True,
    type=click.Choice(["rnadiff", "bamqc", "enrichment"]))
@click.option("--enrichment-taxon", type=click.INT,
    #required=True,
    default=0,
    help="a valid taxon identifiers")
@click.option("--enrichment-kegg-name", type=click.STRING,
    default=None,
    help="a valid KEGG name (automatically filled for 9606 (human) and 10090 (mmusculus)")
@click.option("--enrichment-log2-foldchange-cutoff", type=click.FLOAT,
    default=1,
    show_default=True,
    help="remove events with absolute log2 fold change below this value")
@click.option("--enrichment-padj-cutoff", type=click.FLOAT,
    default=0.05,
    show_default=True,
    help="remove events with pvalue abobe this value default (0.05).")
@click.option("--enrichment-biomart", type=click.STRING,
    default=None,
    help="""you may need a biomart mapping of your identifier for the kegg
pathways analysis. If you do not have this file, you can use 'sequana biomart'
command""")
@click.option("--enrichment-go-only", type=click.BOOL,
    default=False,
    is_flag=True,
    help="""to run only panther db enrichment""")
@click.option("--enrichment-kegg-only", type=click.BOOL,
    default=False,
    is_flag=True,
    help="""to run only kegg patways enrichment""")
@click.option("--enrichment-kegg-pathways-directory", type=click.Path(),
    default=None,
    help="""a place where to find the pathways for each organism""")
@click.option("--enrichment-kegg-background", type=click.INT,
    default=None,
    help="""a background for kegg enrichment. If None, set to number of genes found in KEGG""")
@common_logger
def summary(**kwargs):
    """Create a HTML report for various sequana out

    \b
    * rnadiff: the output of RNADiff pipeline
    * enrichment: the output of RNADiff pipeline
    * bamqc

    Example for the enrichment module:

        sequana summary T1vsT0.complete.xls --module enrichment --enrichment-taxon 10090 
        --enrichment-log2-foldchange-cutoff 2 --enrichment-kegg-only

    The KEGG pathways are loaded and it may take time. Once done, they are saved
    in kegg_pathways/organism and be loaded next time:

        sequana summary T1vsT0.complete.xls --module enrichment --enrichment-taxon 10090 
            --enrichment-log2-foldchange-cutoff 2 --enrichment-kegg-only
            --enrichment-kegg-pathways-directory keff_pathways

    """
    name = kwargs['name']
    module = kwargs['module']

    if module == "bamqc":
        from sequana.modules_report.bamqc import BAMQCModule
        report = BAMQCModule(name, "bamqc.html")
    elif module == "rnadiff":
        from sequana.rnadiff import RNADiffResults
        from sequana.modules_report.rnadiff import RNAdiffModule
        data = RNADiffResults(name)
        report = RNAdiffModule(data)
    elif module == "enrichment":
        from sequana.modules_report.enrichment import Enrichment
        taxon = kwargs['enrichment_taxon']
        if taxon == 0:
            logger.error("You must provide a taxon with --enrichment_taxon")
            return
        keggname = kwargs['enrichment_kegg_name']
        params = {"padj": kwargs['enrichment_padj_cutoff'],
                  "log2_fc": kwargs['enrichment_log2_foldchange_cutoff'],
                  "mapper": kwargs['enrichment_biomart'],
                  "kegg_background": kwargs['enrichment_kegg_background'],
                  "preload_directory": kwargs['enrichment_kegg_pathways_directory'],
                }
        filename = kwargs['enrichment_biomart']
        if filename and os.path.exists(filename) is False:
            logger.error("{} does not exists".format(filename))
            sys.exit(1)
        filename = kwargs['enrichment_kegg_pathways_directory']
        if filename and os.path.exists(filename) is False:
            logger.error("{} does not exists".format(filename))
            sys.exit(1)

        report = Enrichment(name, taxon,
            kegg_organism=keggname, 
            enrichment_params=params,
            go_only=kwargs["enrichment_go_only"],
            kegg_only=kwargs["enrichment_kegg_only"], 
            command=" ".join(['sequana'] + sys.argv[1:]))

@main.command()
@click.option("--mart", default="ENSEMBL_MART_ENSEMBL",
    show_default=True,
    help="A valid mart name")
@click.option("--dataset",  required=True,  
    help="A valid dataset name. e.g. mmusculus_gene_ensembl, hsapiens_gene_ensembl")
@click.option("--attributes",  multiple=True,
    default=["ensembl_gene_id","go_id","entrezgene_id","external_gene_name"],
    show_default=True,
    help="A list of valid attributes to look for in the dataset")
@click.option("--output", default=None,
    help="""by default save results into a CSV file named
    biomart_<dataset>_<YEAR>_<MONTH>_<DAY>.csv""")
@common_logger
def biomart(**kwargs):
    """Retrieve information from biomart and save into CSV file

    This command uses BioMart from BioServices to introspect a MART service
    (--mart) and a specific dataset (default to mmusculus_gene_ensembl). Then,
    for all ensembl IDs, it will fetch the requested attributes (--attributes).
    Finally, it saves the CSV file into an output file (--output). This takes
    about 5-10 minutes to retrieve the data depending on the connection.

    """
    print(kwargs)
    logger.level = kwargs["logger"]

    mart = kwargs['mart']
    attributes = kwargs['attributes']
    dataset = kwargs["dataset"]

    from sequana.enrichment import Mart
    conv = Mart(dataset, mart)
    df = conv.query(attributes)
    conv.save(df, filename=kwargs['output'])


@main.command()
@click.option("-i", "--input", required=True,
    help="The salmon input file.")
@click.option("-o", "--output", required=True,
    help="The feature counts output file")
@click.option("-f", "--gff", required=True,
    help="A GFF file compatible with your salmon file")
@click.option("-a", "--attribute", default="ID",
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


@main.command()
@click.option("-i", "--input", required=True)
@click.option("-o", "--output", required=True)
def gtf_fixer(**kwargs):
    """Reads GTF and fix known issues"""
    from sequana.gtf import GTFFixer
    gtf = GTFFixer(kwargs['input'])
    res = gtf.fix(kwargs['output'])
    print(res)



