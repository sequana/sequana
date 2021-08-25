# -*- coding: utf-8 -*-
#
#  This file is part of Sequana software
#
#  Copyright (c) 2016-2021 - Sequana Development Team
#
#  Distributed under the terms of the 3-clause BSD license.
#  The full license is in the LICENSE file, distributed with this software.
#
#  website: https://github.com/sequana/sequana
#  documentation: http://sequana.readthedocs.io
#
##############################################################################
import functools
import glob
import os
from pathlib import Path
import pkg_resources
import shutil
import subprocess
import sys
import tempfile

import click

#import click_completion
#click_completion.init()

from sequana.utils import config
from sequana import version
from sequana.iem import IEM
from sequana import GFF3
from sequana import FastQ


import colorlog

logger = colorlog.getLogger(__name__)

__all__ = ["main"]


# This can be used by all commands as a simple decorator
def common_logger(func):
    @click.option(
        "--logger",
        default="INFO",
        type=click.Choice(["INFO", "DEBUG", "WARNING", "CRITICAL", "ERROR"]),
    )
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        return func(*args, **kwargs)

    return wrapper


def get_env_vars(ctx, args, incomplete):
    return [k for k in os.environ.keys() if incomplete in k]


CONTEXT_SETTINGS = dict(help_option_names=["-h", "--help"])


pipelines = [
    item.key for item in pkg_resources.working_set if item.key.startswith("sequana")
]
if len(pipelines):
    version += "\nThe following pipelines are installed:\n"
for item in pkg_resources.working_set:
    if item.key.startswith("sequana") and item.key != "sequana":
        version += "\n - {} version: {}".format(item.key, item.version)


@click.group(context_settings=CONTEXT_SETTINGS)
@click.version_option(version=version)
def main(**kwargs):
    """This is the main entry point for a set of Sequana applications.

    Pipelines such as sequana_rnaseq, sequana_variant_calling have their own
    application and help.

    In addition, more advanced tools such as sequana_taxonomy or
    sequana_coverage have their own standalone.


    To setup completion, type this command depending on your shell (bash):

    \b
        eval "$(_SEQUANA_COMPLETE=source_bash sequana)"
        eval "$(_SEQUANA_COMPLETE=source_zsh sequana)"
        eval (env _SEQUANA_COMPLETE=source_fish sequana)

    """
    pass


@main.command()
@click.argument("filename", type=click.STRING, nargs=-1)
@click.option(
    "-o",
    "--output",
    help="filename where to save results. to be used with --head, --tail",
)
@click.option("--count-reads", is_flag=True)
@click.option("--head", type=click.INT, help="number of reads to extract from the head")
@click.option("--merge", is_flag=True)
@click.option("--tail", type=click.INT, help="number of reads to extract from the tail")
def fastq(**kwargs):
    """Set of useful utilities for FastQ manipulation.

    Input file can be gzipped or not. The --output-file

    """

    filenames = kwargs["filename"]
    # users may provide a wildcards such as "A*gz" or list of files.
    if len(filenames) == 1:
        # if existing files or glob, a glob would give the same answer.
        filenames = glob.glob(filenames[0])
    for filename in filenames:
        os.path.exists(filename)

    # could be simplified calling count_reads only once
    if kwargs["count_reads"]:
        for filename in filenames:
            f = FastQ(filename)
            Nreads = f.count_reads()
            Nlines = Nreads * 4
            print(f"Number of reads in {filename}: {Nreads}")
            print(f"Number of lines in {filename}: {Nlines}")
    elif kwargs["head"]:
        for filename in filenames:
            f = FastQ(filename)
            if kwargs["output"] is None:
                logger.error("Please use --output to tell us where to save the results")
                sys.exit(1)
            N = kwargs["head"] * 4
            f.extract_head(N=N, output_filename=kwargs["output"])
    elif kwargs["tail"]:  # pragma: no cover
        raise NotImplementedError
    elif kwargs["merge"]:
        # merge all input files (assuming gz extension)
        extensions = [filename.split(".")[-1] for filename in filenames]
        if set(extensions) != set(["gz"]):
            raise ValueError("Your input FastQ files must be zipped")
        output_filename = kwargs["output"]
        if output_filename is None:
            logger.error("You must use --output filename.gz")
            sys.exit(1)
        if output_filename.endswith(".gz") is False:
            raise ValueError("your output file must end in .gz")

        p1 = subprocess.Popen(["zcat"] + list(filenames), stdout=subprocess.PIPE)
        fout = open(output_filename, "wb")
        p2 = subprocess.run(["pigz"], stdin=p1.stdout, stdout=fout)

    else:  # pragma: no cover
        print("Use one of the commands")


@main.command()
@click.argument("name", type=click.STRING)
@click.option("--check", is_flag=True)
@click.option("--extract-adapters", is_flag=True)
@click.option("--quick-fix", is_flag=True)
@click.option("--output", default=None)
def samplesheet(**kwargs):
    """Utilities to manipulate sample sheet"""
    name = kwargs["name"]
    if kwargs["check"]:
        iem = IEM(name)
        iem.validate()
        logger.info("SampleSheet looks correct")
    elif kwargs["extract_adapters"]:
        iem = IEM(name)
        iem.to_fasta()
    elif kwargs["quick_fix"]:
        iem = IEM(name, tryme=True)
        if kwargs["output"]:
            filename = kwargs["output"]
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
@click.option(
    "--module",
    required=False,
    type=click.Choice(["bamqc", "bam", "fasta", "fastq", "gff"]),
)
def summary(**kwargs):
    """Create a HTML report for various type of NGS formats.

    \b
    * bamqc
    * fastq

    This will process all files in the given pattern (in back quotes)
    sequentially and procude one HTML file per input file.


    Other module all work in the same way. For example, for FastQ files::

        sequana summary one_input.fastq
        sequana summary `ls *fastq`


    """
    names = kwargs["name"]
    module = kwargs["module"]

    if module is None:
        if names[0].endswith("fastq.gz") or names[0].endswith(".fastq"):
            module = "fastq"
        elif names[0].endswith(".bam"):
            module = "bam"
        elif names[0].endswith(".gff") or names[0].endswith("gff3"):
            module = "gff"
        elif names[0].endswith("fasta.gz") or names[0].endswith(".fasta"):
            module = "fasta"
        else:
            logger.error("please use --module to tell us about the input fimes")
            sys.exit(1)

    if module == "bamqc":
        for name in names:
            print(f"Processing {name}")
            from sequana.modules_report.bamqc import BAMQCModule

            report = BAMQCModule(name, "bamqc.html")
    elif (
        module == "fasta"
    ):  # there is no module per se. HEre we just call FastA.summary()
        from sequana.fasta import FastA

        for name in names:
            f = FastA(name)
            f.summary()
    elif (
        module == "fastq"
    ):  # there is no module per se. HEre we just call FastA.summary()
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
    elif module == "gff":
        import pandas as pd
        from sequana import GFF3

        for filename in names:
            ff = GFF3(filename)
            print("#filename: {}".format(filename))
            print("#Number of entries per genetic type:")
            print(ff.df.value_counts("type").to_string())
            print("#Number of duplicated attribute (if any) per attribute:")
            ff.get_duplicated_attributes_per_type()


@main.command()
@click.option(
    "--file1",
    type=click.Path(),
    default=None,
    required=True,
    help="""The first input RNA-seq table to compare""",
)
@click.option(
    "--file2",
    type=click.Path(),
    default=None,
    required=True,
    help="""The second input RNA-seq table to compare""",
)
@common_logger
def rnaseq_compare(**kwargs):
    """Compare 2 tables created by the 'sequana rnadiff' command"""
    from sequana.compare import RNADiffCompare
    from pylab import savefig

    c = RNADiffCompare(kwargs["file1"], kwargs["file2"])
    print(c.r1.summary())
    print(c.r2.summary())
    c.plot_volcano_differences()
    savefig("sequana_rnaseq_compare_volcano.png", dpi=200)


# a dynamic call back function to introspect the design batch column
from sequana.rnadiff import RNADiffAnalysis, RNADesign


def rnadiff_auto_batch_column(ctx, args, incomplete):
    if "--design" in args:
        dfile = args[args.index("--design")]
        d = RNADesign(dfile)
    else:
        d = RNADesign("design.csv")

    batch = (x for x in d.df.columns if x not in {"label", "condition"})
    if len(batch) == 0:
        logger.warning("No batch effect included in your design file")
    else:
        return [c for c in batch if incomplete in c[0]]


@main.command()
@click.option(
    "--annotation",
    type=click.Path(),
    default=None,
    help="""The annotation GFF file used to perform the feature count""",
)
@click.option(
    "--output-directory",
    type=click.Path(),
    default="rnadiff",
    help="""Output directory where are saved the results. Use --force if it exists already""",
)
@click.option(
    "--force/--no-force",
    default=False,
    help="If output directory exists, use this option to erase previous results",
)
@click.option(
    "--features",
    type=click.Path(),
    default="all_features.out",
    help="""The Counts from feature counts. This should be the output of the
            sequana_rnaseq pipeline all_features.out """,
)
# FIXME I think it would be better to have a single file with multiple columns
# for alternative condition (specified using the "condition" option)
@click.option(
    "--design",
    type=click.Path(),
    default="design.csv",
    help="""It should have been generated by sequana_rnaseq. If
not, it must be a comma separated file with two columns. One for the label to be
found in the --features file and one column with the condition to which it
belong. Extra columns can be added to add batch effet. With 3 replicates and 2 conditions,
it should look like:

\b
label,condition
WT1,WT
WT2,WT
WT3,WT
file1,cond1
fileother,cond1
""",
)
@click.option(
    "--condition",
    type=str,
    default="condition",
    help="""The name of the column in design.csv to use as condition
for the differential analysis. Default is 'condition'""",
)
@click.option(
    "--feature-name",
    default="gene",
    help="The feature name compatible with your GFF (default is 'gene')",
)
@click.option(
    "--attribute-name",
    default="ID",
    help="""The attribute used as identifier. Compatible with your GFF (default is 'ID')""",
)
@click.option(
    "--reference",
    type=click.Path(),
    default=None,
    help="""The reference to test DGE against. If provided, conditions not
            involving the reference are ignored. Otherwise all combinations are
            tested""",
)
@click.option(
    "--comparisons",
    type=click.Path(),
    default=None,
    help="""By default, if a reference is provided, all conditions versus that
reference are tested. If no reference, the entire combinatory is performed
(Ncondition * (Ncondition-1) / 2. In both case all condtions found in the
design file are used. If a comparison file is provided, only conditions found in
it will be used. """,
)
@click.option(
    "--cooks-cutoff",
    type=click.Path(),
    default=None,
    help="""if none, let DESeq2 choose the cutoff""",
)
@click.option(
    "--independent-filtering/--no-independent-filtering",
    default=False,
    help="""Do not perform independent_filtering by default. low counts may not
have adjusted pvalues otherwise""",
)
@click.option(
    "--beta-prior/--no-beta-prior",
    default=False,
    help="Use beta prior or not. Default is no beta prior"
)
@click.option(
    "--batch",
    type=str,
    default=None,
    help="""set the column name (in your design) corresponding to the batch
effect to be included in the statistical model as batch ~ condition""",
    autocompletion=rnadiff_auto_batch_column,
)
@click.option(
    "--fit-type",
    default="parametric",
    help="DESeq2 type of fit. Default is 'parametric'",
)
@click.option(
    "--minimum-mean-reads-per-gene",
    default=0,
    help="Filter out fene where the mean number of reads is below this value. By default all genes are kept",
)
@click.option(
    "--keep-all-conditions/--no-keep-all-conditions",
    default=False,
    help="""Even though sub set of comparisons are provided, keep all conditions
in the analysis and report only the provided comparisons"""
)
@click.option(
    "--hover-name",
    default=None,
    help="""In volcano plot, we set the hover name to Name if present in the GFF,
otherwise to gene_id if present, then locus_tag, and finally ID and gene_name. One can specify
a hover name to be used with this option""",
)
@common_logger
def rnadiff(**kwargs):
    """Perform RNA-seq differential analysis and reporting.

        This command performs the differential analysis of gene expression. The
        analysis is performed on feature counts generated by a RNA-seq analysis
        (see e.g. https://github.com/sequana/rnaseq pipeline). The analysis is
        performed by DESeq2. A HTML report is created as well as a set of output
        files, including summary table of the analysis.

        To perform this analysis, you will need the GFF file used during the RNA-seq
        analysis, the feature stored altogether in a single file, an experimental
        design file, and the feature and attribute used during the feature count.

        Here is an example:

    \b
            sequana rnadiff --annotation Lepto.gff
                --design design.csv --features all_features.out
                 --feature-name gene --attribute-name ID

        Batch effet can be included by adding a column in the design.csv file. For
        example if called 'day', you can take this information into account using
        '--batch day'


    """
    import pandas as pd
    from sequana.featurecounts import FeatureCount
    from sequana.rnadiff import RNADiffAnalysis, RNADesign
    from sequana.modules_report.rnadiff import RNAdiffModule

    # FIXME if we use the colorlog logger here, it does not work....
    # FIXME not clear why so we import again the sequana logger
    # FIXME looks like it works everywhere but not in the main script
    from sequana import logger
    # seems to be fixed in aug 2021

    logger.setLevel(kwargs["logger"])

    from easydev import cmd_exists, mkdirs
    if not cmd_exists("Rscript"):
        logger.critical("""Rscript not found; You will need R and the DESeq2 package to be installed. 
You may install it yourself or use damona using the rtools:1.0.0 image """)
        sys.exit(1)

    outdir = kwargs["output_directory"]
    feature = kwargs["feature_name"]
    attribute = kwargs["attribute_name"]

    if os.path.exists(outdir) and not kwargs["force"]:
        logger.error(f"{outdir} exist already. Use --force to overwrite")
        sys.exit(1)

    if kwargs["annotation"]:
        gff_filename = kwargs["annotation"]
        logger.info(f"Checking annotation file (feature and attribute)")
        gff = GFF3(gff_filename)
        if feature not in gff.features:
            logger.error(
                f"{feature} not found in the GFF. Most probably a wrong feature name"
            )
            sys.exit(1)
        attributes = gff.get_attributes(feature)
        if attribute not in attributes:
            logger.error(
                f"{attribute} not found in the GFF for the provided feature. Most probably a wrong feature name. Please change --attribute-name option or do not provide any GFF"
            )
            sys.exit(1)
    else:
        gff = None

    comparisons = kwargs["comparisons"]
    if comparisons:
        compa_df = pd.read_csv(comparisons)
        comparisons = list(zip(compa_df["alternative"], compa_df["reference"]))

    logger.info(f"Differential analysis to be saved into ./{outdir}")
    for k in sorted(
        [
            "independent_filtering",
            "beta_prior",
            "batch",
            "cooks_cutoff",
            "fit_type",
            "reference",
        ]
    ):
        logger.info(f"  Parameter {k} set to : {kwargs[k]}")

    # The readl analysis is here
    r = RNADiffAnalysis(
        kwargs["features"],
        kwargs["design"],
        kwargs["condition"],
        keep_all_conditions=kwargs['keep_all_conditions'],
        batch=kwargs["batch"],
        comparisons=comparisons,
        reference=kwargs['reference'],
        fc_feature=feature,
        fc_attribute=attribute,
        outdir=outdir,
        gff=gff,
        cooks_cutoff=kwargs.get("cooks_cutoff"),
        independent_filtering=kwargs.get("independent_filtering"),
        beta_prior=kwargs.get("beta_prior"),
        fit_type=kwargs.get("fit_type"),
        minimum_mean_reads_per_gene=kwargs.get("minimum_mean_reads_per_gene"),
    )

    try:
        logger.info(f"Running DGE. Saving results into {outdir}")
        results = r.run()
        results.to_csv(f"{outdir}/rnadiff.csv")
    except Exception as err:
        logger.error(err)
        logger.error(f"please see {outdir}/code/rnadiff.err file for errors")
        sys.exit(1)

    logger.info(f"Reporting. Saving in summary.html")

    # this define the output directory where summary.html is saved
    config.output_dir = outdir

    import seaborn as sns

    report = RNAdiffModule(
        outdir,
        gff=gff,
        fc_attribute=attribute,
        fc_feature=feature,
        alpha=0.05,
        log2_fc=0,
        condition=kwargs["condition"],
        annot_cols=None,
        pattern="*vs*_degs_DESeq2.csv",
        palette=sns.color_palette(desat=0.6, n_colors=13),
        hover_name=kwargs["hover_name"],
    )

    #
    # save info.txt with sequana version
    teardown(outdir)


@main.command()
@click.option(
    "--mart",
    default="ENSEMBL_MART_ENSEMBL",
    show_default=True,
    help="A valid mart name",
)
@click.option(
    "--dataset",
    required=True,
    help="A valid dataset name. e.g. mmusculus_gene_ensembl, hsapiens_gene_ensembl",
)
@click.option(
    "--attributes",
    default="ensembl_gene_id,go_id,entrezgene_id,external_gene_name",
    show_default=True,
    help="""A valid set of attributes to look for in the dataset. Multiple
attributes are separeted by a comma (no spaces accepted)""",
)
@click.option(
    "--output",
    default=None,
    help="""by default save results into a CSV file named
    biomart_<dataset>_<YEAR>_<MONTH>_<DAY>.csv""",
)
@common_logger
def biomart(**kwargs):
    """Retrieve information from biomart and save into CSV file

    This command uses BioMart from BioServices to introspect a MART service
    (--mart) and a specific dataset (default to mmusculus_gene_ensembl). Then,
    for all ensembl IDs, it will fetch the requested attributes (--attributes).
    Finally, it saves the CSV file into an output file (--output). This takes
    about 5-10 minutes to retrieve the data depending on the connection.

    Example:

        sequana biomart --mart mmusculus_gene_ensembl --mart ENSEMBL_MART_ENSEMBL \
            --dataset mmusculus_gene_ensembl \
            --attributes ensembl_gene_id,external_gene_name,go_id \
            --output test.csv

    """
    logger.setLevel(kwargs["logger"])

    mart = kwargs["mart"]
    attributes = kwargs["attributes"]
    dataset = kwargs["dataset"]

    from sequana.enrichment.mart import Mart

    conv = Mart(dataset, mart)
    df = conv.query(attributes.split(","))
    conv.save(df, filename=kwargs["output"])


@main.command()
@click.option("-i", "--input", required=True, help="The salmon input file.")
@click.option("-o", "--output", required=True, help="The feature counts output file")
@click.option(
    "-f", "--gff", required=True, help="A GFF file compatible with your salmon file"
)
@click.option(
    "-a",
    "--attribute",
    default="ID",
    help="A valid attribute to be found in the GFF file and salmon input",
)
@click.option("-a", "--feature", default="gene", help="A valid feature")
def salmon(**kwargs):
    """Convert output of Salmon into a feature counts file """
    from sequana import salmon

    salmon_input = kwargs["input"]
    output = kwargs["output"]
    if os.path.exists(salmon_input) is False:
        logger.critical("Input file does not exists ({})".format(salmon_input))
    gff = kwargs["gff"]
    attribute = kwargs["attribute"]
    feature = kwargs["feature"]

    # reads file generated by salmon and generated count file as expected by
    # DGE.
    s = salmon.Salmon(salmon_input, gff)
    s.save_feature_counts(output, feature=feature, attribute=attribute)


@main.command()
@click.option("-i", "--input", required=True)
@click.option("-o", "--output", required=True)
def gtf_fixer(**kwargs):
    """Reads GTF and fix known issues (exon and genes uniqueness)"""
    from sequana.gtf import GTFFixer

    gtf = GTFFixer(kwargs["input"])
    res = gtf.fix_exons_uniqueness(kwargs["output"])
    # res = gtf.fix_exons_uniqueness(kwargs['output'])
    print(res)


# This will be a complex command to provide HTML summary page for
# input files (e.g. bam), or results from pipelines. For each module,
# we should have corresponding option that starts with the module's name
# This can also takes as input various types of data (e.g. FastA)
@main.command()
@click.argument("name", type=click.Path(exists=True), nargs=1)
@click.option(
    "--annotation-attribute",
    type=click.STRING,
    # required=True,
    default="Name",
    help="a valid taxon identifiers",
)
@click.option(
    "--panther-taxon",
    type=click.INT,
    # required=True,
    default=0,
    help="a valid taxon identifiers",
)
@click.option(
    "--kegg-name",
    type=click.STRING,
    default=None,
    help="a valid KEGG name (automatically filled for 9606 (human) and 10090 (mmusculus)",
)
@click.option(
    "--log2-foldchange-cutoff",
    type=click.FLOAT,
    default=1,
    show_default=True,
    help="remove events with absolute log2 fold change below this value",
)
@click.option(
    "--padj-cutoff",
    type=click.FLOAT,
    default=0.05,
    show_default=True,
    help="remove events with pvalue abobe this value default (0.05).",
)
@click.option(
    "--biomart",
    type=click.STRING,
    default=None,
    help="""you may need a biomart mapping of your identifier for the kegg
pathways analysis. If you do not have this file, you can use 'sequana biomart'
command""",
)
@click.option(
    "--go-only",
    type=click.BOOL,
    default=False,
    is_flag=True,
    help="""to run only panther db enrichment""",
)
@click.option(
    "--plot-linearx",
    type=click.BOOL,
    default=False,
    is_flag=True,
    help="""Default is log2 fold enrichment in the plots. use this to use linear scale""",
)
@click.option(
    "--compute-levels",
    type=click.BOOL,
    default=False,
    is_flag=True,
    help="""to compute the GO levels (slow) in the plots""",
)
@click.option(
    "--max-genes",
    type=click.INT,
    default=2000,
    help="""Maximum number of genes (up or down) to use in PantherDB, which is limited to about 3000""",
)
@click.option(
    "--kegg-only",
    type=click.BOOL,
    default=False,
    is_flag=True,
    help="""to run only kegg patways enrichment""",
)
@click.option(
    "--kegg-pathways-directory",
    type=click.Path(),
    default=None,
    help="""a place where to find the pathways for each organism""",
)
@click.option(
    "--kegg-background",
    type=click.INT,
    default=None,
    help="""a background for kegg enrichment. If None, set to number of genes found in KEGG""",
)
@common_logger
def enrichment(**kwargs):
    """Create a HTML report for various sequana out

    \b
    * enrichment: the output of RNADiff pipeline

    Example for the enrichment module:

        sequana enrichment rnadiff.csv --panther-taxon 10090
            --log2-foldchange-cutoff 2 --kegg-only

    The KEGG pathways are loaded and it may take time. Once done, they are saved
    in kegg_pathways/organism and be loaded next time:

        sequana enrichment rnadiff/rnadiff.csv
            --panther-taxon 189518 \
            --log2-foldchange-cutoff 2 --kegg-only \
            --kegg-name lbi\
            --annotation file.gff

    """
    import pandas as pd
    from sequana.modules_report.enrichment import Enrichment

    logger.setLevel(kwargs["logger"])

    taxon = kwargs["panther_taxon"]
    if taxon == 0:
        logger.error("You must provide a taxon with --panther-taxon")
        sys.exit(1)
    keggname = kwargs["kegg_name"]
    params = {
        "padj": kwargs["padj_cutoff"],
        "log2_fc": kwargs["log2_foldchange_cutoff"],
        "max_entries": kwargs["max_genes"],
        "mapper": kwargs["biomart"],
        "kegg_background": kwargs["kegg_background"],
        "preload_directory": kwargs["kegg_pathways_directory"],
        "plot_logx": not kwargs["plot_linearx"],
        "plot_compute_levels": kwargs["compute_levels"],
    }
    filename = kwargs["biomart"]
    if filename and os.path.exists(filename) is False:
        logger.error("{} does not exists".format(filename))
        sys.exit(1)
    filename = kwargs["kegg_pathways_directory"]
    if filename and os.path.exists(filename) is False:
        logger.error("{} does not exists".format(filename))
        sys.exit(1)

    rnadiff_file = kwargs["name"]
    logger.info(f"Reading {rnadiff_file}")
    rnadiff = pd.read_csv(rnadiff_file, index_col=0, header=[0, 1])

    # now that we have loaded all results from a rnadiff analysis, let us
    # perform the enrichment for each comparison found in the file
    annot_col = kwargs["annotation_attribute"]
    Nmax = kwargs["max_genes"]

    from sequana.utils import config

    for compa in rnadiff.columns.levels[0]:
        if compa not in ["statistics", "annotation"]:
            # get gene list
            df = rnadiff[compa].copy()

            # we add the annotation
            for x in rnadiff["annotation"].columns:
                df[x] = rnadiff["annotation"][x]

            # now we find the gene lists
            padj = params["padj"]
            log2fc = params["log2_fc"]
            df = df.query(
                "(log2FoldChange >=@log2fc or log2FoldChange<=-@log2fc) and padj <= @padj"
            )
            df.reset_index(inplace=True)

            dfup = df.sort_values("log2FoldChange", ascending=False)
            up_genes = list(dfup.query("log2FoldChange > 0")[annot_col].fillna("undefined"))[:Nmax]

            dfdown = df.sort_values("log2FoldChange", ascending=True)
            down_genes = list(dfdown.query("log2FoldChange < 0")[annot_col].fillna("undefined"))[:Nmax]

            all_genes = list(
                df.sort_values("log2FoldChange", key=abs, ascending=False)[annot_col].dropna())[:Nmax]

            gene_dict = {
                "up": up_genes,
                "down": down_genes,
                "all": all_genes,
            }

            Nup = len(up_genes)
            Ndown = len(down_genes)
            N = Nup + Ndown
            logger.info(f"Computing enrichment for the {compa} case")
            logger.info(
                f"Found {Nup} genes up-regulated, {Ndown} down regulated ({N} in total) based on padj<{padj} and absolute log2 fold change < {log2fc}. Note that if the attribute is not set, some genes may be ignored"
            )
            config.output_dir = f"enrichment/{compa}"

            if not os.path.exists("enrichment"):
                os.mkdir("enrichment")

            Enrichment(
                gene_dict,
                taxon,
                df,
                kegg_organism=keggname,
                enrichment_params=params,
                go_only=kwargs["go_only"],
                kegg_only=kwargs["kegg_only"],
                command=" ".join(["sequana"] + sys.argv[1:]),
            )


@main.command()
@click.option(
    "--search-kegg",
    type=click.Path(),
    default=None,
    help="""Search a pattern amongst all KEGG organisms""",
)
@click.option(
    "--search-panther",
    type=click.Path(),
    default=None,
    help="""Search a pattern amongst all KEGG organism""",
)
@common_logger
def taxonomy(**kwargs):
    """Tool to retrieve taxonomic information.

    sequana taxonomy --search-kegg leptospira
    """

    if kwargs["search_kegg"]:
        from sequana.kegg import KEGGHelper

        k = KEGGHelper()
        results = k.search(kwargs["search_kegg"].lower())
        print(results)
    elif kwargs["search_panther"]:
        import pandas as pd
        from sequana import sequana_data

        df = pd.read_csv(sequana_data("panther.csv"), index_col=0)

        pattern = kwargs["search_panther"]
        f1 = df[[True if pattern in x else False for x in df["name"]]]
        f2 = df[[True if pattern in x else False for x in df.short_name]]
        f3 = df[[True if pattern in x else False for x in df.long_name]]
        indices = list(f1.index) + list(f2.index) + list(f3.index)

        if len(indices) == 0:
            # maybe it is a taxon ID ?
            f4 = df[[True if pattern in str(x) else False for x in df.taxon_id]]
            indices = list(f4.index)
        indices = set(indices)
        print(df.loc[indices])

@main.command()
@click.argument("input", type=click.Path(exists=True))
@click.argument("output")
@click.option(
    "--features",
    type=click.Path(),
    default="gene",
    help="""list of features to be  extracted""",
)
@common_logger
def gff_to_light_gff(**kwargs):
    """Extract the feature of interest in the input GFF to create a light version

    sequana gff-to-light-gff input.gff output.gff --features gene,exon 

    """
    from sequana import logger
    logger.setLevel(kwargs["logger"])
    filename = kwargs["input"]
    assert filename.endswith(".gff") or filename.endswith(".gff3")

    g = GFF3(filename)
    g.read_and_save_selected_features(kwargs["output"], features=kwargs["features"])


@main.command()
@click.argument("gff_filename", type=click.Path(exists=True))
@common_logger
def gff_to_gtf(**kwargs):
    """Convert a GFF file into GTF

    This is experimental convertion. Use with care.

    """
    filename = kwargs["gff_filename"]
    assert filename.endswith(".gff") or filename.endswith(".gff3")

    g = GFF3(filename)
    if filename.endswith(".gff"):
        g.to_gtf(os.path.basename(filename).replace(".gff", ".gtf"))
    elif filename.endswith(".gff3"):
        g.to_gtf(os.path.basename(filename).replace(".gff3", ".gtf"))


def teardown(workdir):
    # common function to be used by subcommands to store called command
    from pathlib import Path
    from easydev import mkdirs

    workdir = Path(workdir)
    mkdirs(workdir / ".sequana")
    with open(Path(workdir) / ".sequana" / "info.txt", "w") as fout:
        from sequana import version

        fout.write(f"# sequana version: {version}\n")
        fout.write(" ".join(["sequana"] + sys.argv[1:]))
