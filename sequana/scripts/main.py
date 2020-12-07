#-*- coding: utf-8 -*-

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
        version += "\n - {} version: {}".format(item.key, item.version)


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
@click.argument('filename', type=click.STRING, nargs=-1)
@click.option("-o", "--output",
    help="filename where to save results. to be used with --head, --tail")

@click.option("--count-reads", is_flag=True)
@click.option("--head", type=click.INT,
    help='number of reads to extract from the head')
@click.option("--merge", is_flag=True)
@click.option("--tail", type=click.INT,
    help="number of reads to extract from the tail")
def fastq(**kwargs):
    """Set of useful utilities for FastQ manipulation.

    Input file can be gzipped or not. The --output-file

    """
    from sequana.fastq import FastQ

    filenames = kwargs['filename']
    # users may provide a wildcards such as "A*gz" or list of files.
    if len(filenames) == 1:
        # if existing files or glob, a glob would give the same answer.
        filenames = glob.glob(filenames[0])
    for filename in filenames:
        os.path.exists(filename)

    # could be simplified calling count_reads only once
    if kwargs['count_reads']:
        for filename in filenames:
            f = FastQ(filename)
            Nreads = f.count_reads()
            Nlines = Nreads * 4
            print(f"Number of reads in {filename}: {Nreads}")
            print(f"Number of lines in {filename}: {Nlines}")
    elif kwargs['head']:
        for filename in filenames:
            f = FastQ(filename)
            if kwargs['output'] is None:
                logger.error("Please use --output to tell us where to save the results")
                sys.exit(1)
            N = kwargs['head'] * 4
            f.extract_head(N=N, output_filename=kwargs['output'])
    elif kwargs['tail']: #pragma: no cover
        raise NotImplementedError
    elif kwargs['merge']:
        import subprocess
        # merge all input files (assuming gz extension)
        extensions = [filename.split(".")[-1] for filename in filenames]
        if set(extensions) != set(['gz']):
            raise ValueError("Your input FastQ files must be zipped")
        output_filename = kwargs['output']
        if output_filename is None:
            logger.error("You must use --output filename.gz")
            sys.exit(1)
        if output_filename.endswith(".gz") is False:
            raise ValueError("your output file must end in .gz")

        p1 = subprocess.Popen(['zcat'] +  list(filenames),  stdout=subprocess.PIPE)
        fout = open(output_filename, 'wb')
        p2 = subprocess.run(['pigz'], stdin=p1.stdout, stdout=fout)

    else:  #pragma: no cover
        print("Use one of the commands")


@main.command()
@click.argument('name', type=click.STRING)
@click.option('--check', is_flag=True)
@click.option('--extract-adapters', is_flag=True)
@click.option('--quick-fix', is_flag=True)
@click.option('--output', default=None)
def samplesheet(**kwargs):
    """Utilities to manipulate sample sheet"""
    name = kwargs['name']
    from sequana.iem import IEM
    if kwargs['check']:
        iem = IEM(name)
        iem.validate()
        logger.info("SampleSheet looks correct")
    elif kwargs["extract_adapters"]:
        iem = IEM(name)
        iem.to_fasta()
    elif kwargs["quick_fix"]:
        iem = IEM(name, tryme=True)
        if kwargs['output']:
            filename = kwargs['output']
        else:
            filename = name + ".fixed"
        logger.info("Saving fixed version in {}".format(filename))
        iem.quick_fix(output_filename=filename)


# This will be a complex command to provide HTML summary page for
# input files (e.g. bam), or results from pipelines. For each module,
# we should have corresponding option that starts with the module's name
# This can also takes as input various types of data (e.g. FastA)
@main.command()
@click.argument("name", type=click.Path(exists=True), nargs=-1)
@click.option("--module",
    required=False,
    type=click.Choice(["rnadiff", "bamqc", "bam", "enrichment", "fasta", "fastq"]))
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
@click.option("--enrichment-plot-linearx", type=click.BOOL,
    default=False,
    is_flag=True,
    help="""Default is log2 fold enrichment in the plots. use this to use linear scale""")
@click.option("--enrichment-compute-levels", type=click.BOOL,
    default=False,
    is_flag=True,
    help="""to compute the GO levels (slow) in the plots""")
@click.option("--enrichment-max-genes", type=click.INT,
    default=3000,
    help="""Maximum number of genes (up or down) to use in PantherDB, which is limited to about 3000""")
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
    * fastq

    Example for the enrichment module:

        sequana summary T1vsT0.complete.xls --module enrichment --enrichment-taxon 10090 
        --enrichment-log2-foldchange-cutoff 2 --enrichment-kegg-only

    The KEGG pathways are loaded and it may take time. Once done, they are saved
    in kegg_pathways/organism and be loaded next time:

        sequana summary T1vsT0.complete.xls --module enrichment --enrichment-taxon 10090 
            --enrichment-log2-foldchange-cutoff 2 --enrichment-kegg-only
            --enrichment-kegg-pathways-directory kegg_pathways

    """
    names = kwargs['name']
    module = kwargs['module']

    if module is None:
        if names[0].endswith('fastq.gz') or names[0].endswith('.fastq'):
            module = "fastq"
        elif  names[0].endswith('.bam'):
            module = "bam"
        elif names[0].endswith('fasta.gz') or names[0].endswith('.fasta'):
            module = "fasta"
        else:
            logger.error("please use --module to tell us about the input fimes")
            sys.exit(1)

    if module == "bamqc":
        for name in names:
            print(f"Processing {name}")
            from sequana.modules_report.bamqc import BAMQCModule
            report = BAMQCModule(name, "bamqc.html")
    elif module == "rnadiff":
        for name in names:
            from sequana.rnadiff import RNADiffResults
            from sequana.modules_report.rnadiff import RNAdiffModule
            if name.split(".")[-1] == 'xls':
                output_filename = name.replace(".xls", ".html")
            else:
                output_filename = name + ".html"
            output_filename = output_filename.split("/")[-1]
            logger.info(f"Processing {name} into {output_filename}")
            data = RNADiffResults(name)
            report = RNAdiffModule(data, output_filename)
    elif module == "enrichment":
        try:
            name = names[0]
        except:
            logger.error()
            sys.exit(1)
        from sequana.modules_report.enrichment import Enrichment
        taxon = kwargs['enrichment_taxon']
        if taxon == 0:
            logger.error("You must provide a taxon with --enrichment_taxon")
            return
        keggname = kwargs['enrichment_kegg_name']
        params = {"padj": kwargs['enrichment_padj_cutoff'],
                  "log2_fc": kwargs['enrichment_log2_foldchange_cutoff'],
                  "max_entries": kwargs['enrichment_max_genes'],
                  "mapper": kwargs['enrichment_biomart'],
                  "kegg_background": kwargs['enrichment_kegg_background'],
                  "preload_directory": kwargs['enrichment_kegg_pathways_directory'],
                  "plot_logx": not kwargs['enrichment_plot_linearx'],
                  "plot_compute_levels": kwargs['enrichment_compute_levels'],
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
    elif module == "fasta": # there is no module per se. HEre we just call FastA.summary()
        from sequana.fasta import FastA
        for name in names:
            f = FastA(name)
            f.summary()
    elif module == "fastq": # there is no module per se. HEre we just call FastA.summary()
        from sequana.fastq import FastQ
        from sequana import FastQC
        for filename in names:
            ff = FastQC(filename, max_sample=1e6, verbose=False)
            stats = ff.get_stats()
            print(stats)
    elif module == "bam": 
        import pandas as pd
        from sequana import BAM
        for filename in names:
            ff = BAM(filename)
            stats = ff.get_stats()
            df = pd.Series(stats).to_frame().T
            print(df)




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
    """Reads GTF and fix known issues (exon and genes uniqueness)"""
    from sequana.gtf import GTFFixer
    gtf = GTFFixer(kwargs['input'])
    res = gtf.fix_exons_uniqueness(kwargs['output'])
    #res = gtf.fix_exons_uniqueness(kwargs['output'])
    print(res)



