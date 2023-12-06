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
""".. rubric:: Standalone application dedicated to coverage"""
import argparse
import gc as garbage
import glob
import os
import subprocess
import sys

import colorlog
import rich_click as click
from easydev import AttrDict, shellcmd

from sequana import sequana_data
from sequana import version as sequana_version
from sequana.bedtools import ChromosomeCov, SequanaCoverage
from sequana.modules_report.coverage import ChromosomeCoverageModule, CoverageModule
from sequana.scripts.common import teardown
from sequana.utils import config

from .utils import CONTEXT_SETTINGS

logger = colorlog.getLogger(__name__)


CONTEXT_SETTINGS = dict(help_option_names=["-h", "--help"])


click.rich_click.OPTION_GROUPS = {
    "sequana_coverage": [
        {
            "name": "Required",
            "options": ["--input-file"],
        },
        {
            "name": "Detection Algorithm Options",
            "options": [
                "--circular",
                "--mixture-models",
                "--window-median",
                "--window-gc",
                "--low-threshold",
                "--high-threshold",
                "--threshold",
                "--cnv-clustering",
                "--clustering-parameter",
            ],
        },
        {
            "name": "Modifiers",
            "options": ["--annotation-file", "--reference-file", "--chromosome", "--chunk-size", "--binning"],
        },
        {
            "name": "Download utilities",
            "options": ["--download-reference", "--download-genbank", "--database"],
        },
        {
            "name": "Output files",
            "options": ["--no-multiqc", "--output-directory"],
        },
        {
            "name": "Behaviour",
            "options": ["--version", "--level", "--debug-level", "--help"],
        },
    ],
}


def download_reference(ctx, param, value):
    if value:
        database = ctx.params.get("database", "ENA")
        click.echo(f"Downloading reference {value} from {database}...", nl=False)
        from bioservices.apps import download_fasta as df

        df.download_fasta(value, method="ENA")
        click.echo("done")
        sys.exit(0)
    return value


def download_genbank(ctx, param, value):
    if value:
        click.echo(f"Downloading genbank {value}...", nl=False)
        from sequana.snpeff import download_fasta_and_genbank

        download_fasta_and_genbank(value, value, genbank=True, fasta=False)
        click.echo("done")
        sys.exit(0)
    return value


@click.command(context_settings=CONTEXT_SETTINGS)
@click.option(
    "--database",
    "database",
    default="ENA",
    type=click.Choice(["ENA", "EUtils"]),
    is_eager=True,
    # expose_value=False,
    help="Download the reference from one of these database (default ENA)",
)
@click.option(
    "--download-reference",
    "download_reference",
    default=None,
    type=click.STRING,
    callback=download_reference,
    is_eager=True,
    # expose_value=False,
    help="Provide accession number to download a reference from NCBI/ENA",
)
@click.option(
    "--download-genbank",
    "download_genbank",
    default=None,
    type=click.STRING,
    callback=download_genbank,
    is_eager=True,
    help="Provide accession number to download annotation from NCBI/ENA",
)
@click.option(
    "-i",
    "--input-file",
    "input",
    type=click.Path(),
    required=True,
    help="""Input file in BED or BAM format. If a BAM file is provided, it will be converted locally to a BED file using genomecov, which must be installed.""",
)
@click.option(
    "-a",
    "--annotation-file",
    "annotation",
    type=click.STRING,
    default=None,
    help="a valid gff3 or genbank annotation. Recommended for prokaryotes or small eukaryotes only.",
)
@click.version_option(sequana_version)
@click.option(
    "-o",
    "--circular",
    is_flag=True,
    show_default=True,
    # help="""If the DNA of the organism is circular (typically
    #        viruses or bacteria), set to True"""
)
@click.option(
    "-c",
    "--chromosome",
    "chromosome",
    type=click.STRING,
    default="",
    help="By default all chromosomes are analysed. You may want to analyse only one" " by using this parameter",
)
@click.option(
    "-k",
    "--mixture-models",
    "k",
    type=click.INT,
    help="Number of mixture models to use (default 2, although if sequencing depth is below 8, k is set to 1 automatically). To ignore that behaviour set k to the required value",
    default=2,
)
@click.option(
    "--debug-level",
    "logging_level",
    default="INFO",
    type=click.Choice(["DEBUG", "INFO", "WARNING", "CRITICAL", "ERROR"]),
    help="set to DEBUG, INFO, WARNING, CRITICAL, ERROR",
)
@click.option(
    "--level",
    "logging_level",
    default="INFO",
    type=click.Choice(["DEBUG", "INFO", "WARNING", "CRITICAL", "ERROR"]),
    help="set to DEBUG, INFO, WARNING, CRITICAL, ERROR",
)
@click.option(
    "-g",
    "--window-gc",
    "w_gc",
    type=click.INT,
    default=101,
    help="Length of the running window to compute the GC content",
)
@click.option(
    "--no-multiqc",
    "skip_multiqc",
    is_flag=True,
    show_default=True,
    help="Do not create any multiqc HTML page.",
)
@click.option("--output-directory", "output_directory", default="report", help="name of the output (report) directory.")
@click.option(
    "-r",
    "--reference-file",
    "reference",
    type=click.STRING,
    default=None,
    help="If available, you can provide a reference (ENA/NCBI). It must have the same length as the one used to create the BAM or BED file. If provided, it is used to create the coverage versus GC content image. Recommended for prokaryotes or small eukaryotes only.",
)
@click.option(
    "-w",
    "--window-median",
    "w_median",
    type=click.IntRange(200),
    help="Length of the running median window (default 20,001, recommended for bacteria).  For short genome (below 100000 bases), we set this parameter to one fifth of the genome length. minimal value is 125 (one fourth of 500bp that should be minimal length of contigs)",
    default=20001,
)
@click.option(
    "-L",
    "--low-threshold",
    "low_threshold",
    default=-4,
    type=click.FLOAT,
    help=("lower threshold (zscore) of the confidence interval. " "Overwrite value given by --threshold/-T"),
)
@click.option(
    "-H",
    "--high-threshold",
    "high_threshold",
    default=4,
    type=click.FLOAT,
    help=("higher threshold (zscore) of the confidence interval. " "Overwrite value given by --threshold/-T"),
)
@click.option(
    "-T",
    "--threshold",
    "threshold",
    default=4,
    type=click.FLOAT,
    help="set lower and higher thresholds of the confidence interval. ",
)
@click.option(
    "-C",
    "--clustering-parameter",
    "double_threshold",
    default=0.5,
    type=click.FLOAT,
    help="set lower and higher double threshold parameter (in [0,1]). Do not use value close to zero. Ideally, around 0.5. lower value will tend to cluster more than higher value",
)
@click.option(
    "-s",
    "--chunk-size",
    "chunksize",
    type=click.IntRange(1000000),  # 1e6
    default=10000000,  # 1e7
    help="Length of the chunk to be used for the analysis. ",
)
@click.option(
    "-B",
    "--binning",
    "binning",
    type=click.IntRange(2),
    default=None,
    help="merge consecutive (non overlapping) data points, taking the mean. This is useful for large genome (e.g. human). This allows a faster computation, especially for CNV detection were only large windows are of interest. For instance, using a binning of 50 or 100 allows the human genome to be analysed.",
)
@click.option(
    "--cnv-clustering",
    "cnv_clustering",
    default=-1,
    type=click.INT,
    help="Two consecutive ROIs are merged when their distance in bases is below this parameter. If set to -1, not used. ",
)
def main(**kwargs):
    """Welcome to SEQUANA -- Coverage standalone

        ----

       Extract and plot coverage of one or more chromosomes/contigs from a BED or BAM
       file. In addition, the running median used in conjunction with double thresholds
       extract regions of interests (ROI) for low or high coverage. A reference may be
       provided to plot the coverage versus GC content. An annotation file may be provided
       to annotate the found ROIs.

       The input file should be one of the following:

       - a BED file that is a tabulated file at least 3 columns.
         The first column being the reference, the second is the position
         and the third column contains the coverage itself.
       - or a BAM file that is converted automatically
         into a BED file using the following command:

            bedtools genomecov -d -ibam FILE.bam

       Note that this is different from another command that could perform the conversion:

           samtools depth -aa input.bam > output.bed

       If the reference is provided, an additional plot showing the coverage versus
       GC content is also shown.

       Here are some examples

           sequana_coverage --input-file file.bed --window-median 1001
           sequana_coverage --input-file file.bam --window-median 1001 -r <REFERENCE.fa>

       An other interesting option is to provide a BED file with 4 columns. The
       fourth column being another coverage data created with a filter. One can
       create such a file only from the BAM file using samtools as follows given
       the original unfiltered BAM file named input.bam:

           samtools view -q 35  -o data.filtered.bam input.bam
           samtools depth input.bam data.filtered.bam  -aa > test.bed
           sequana_coverage --input test.bed --show-html

       Note that the first file is the filtered one, and the second file is the
       unfiltered one.

       Note for multi chromosome and genbank features: for now, you will need to call
       sequana_coverage for each chromosome individually since we accept only one
       genbank as input parameter:

           sequana_coverage --input-file file.bed --genbank chrom1.gbk -c 1

       Large genomes:
       --------------

       If your input data is large and does not fit into memory, use the --binning BIN
       options to average data into bin of BIN values.

       CNV cases:
       --------------

       By default, sequana_coverage identify events as small as 1 bin. For the CNV
       detection case, you may want to cluster events. the --cnv-merging DELTA option
       merges consecutives events whose distance is smaller that DELTA


    ----

    AUTHORS: Thomas Cokelaer, Dimitri Desvillechabrol
    Documentation: http://sequana.readthedocs.io
    Issues: http://github.com/sequana/sequana
    """
    options = AttrDict(**kwargs)

    logger.setLevel(options.logging_level)

    if options.annotation:
        assert os.path.exists(options.annotation), "%s does not exists" % options.annotation

    logger.info(f"Reading {options.input}. This may take time ")

    # Convert BAM to BED
    if options.input.endswith(".bam"):
        bedfile = options.input.replace(".bam", ".bed")
        logger.info("Converting BAM into BED file")
        shellcmd(f"bedtools genomecov -d -ibam {options.input} > {bedfile}")
    elif options.input.endswith(".bed"):
        bedfile = options.input
    else:
        raise ValueError("Input file must be a BAM or BED file")

    # Set the thresholds
    if options.low_threshold is None:
        options.low_threshold = -options.threshold

    if options.high_threshold is None:
        options.high_threshold = options.threshold

    # Now we can create the instance of GenomeCoverage
    if "," in options.chromosome:
        chrom_list = options.chromosome.split(",")
    elif options.chromosome:
        chrom_list = [options.chromosome]
    else:
        chrom_list = []

    # initialisation (reading BED to get positions of chromosomes and chromosome names
    gc = SequanaCoverage(
        bedfile,
        options.annotation,
        options.low_threshold,
        options.high_threshold,
        options.double_threshold,
        options.double_threshold,
        chunksize=options.chunksize,
        chromosome_list=chrom_list,
        reference_file=options.reference,
        gc_window_size=options.w_gc,
    )

    # some information fo end-users,
    logger.info("There are %s chromosomes/contigs." % len(gc))
    for this in gc.chrom_names:
        end = gc.positions[this]["end"]
        start = gc.positions[this]["start"]
        data = (this, gc.positions[this]["start"], gc.positions[this]["end"], end - start)
        logger.info("    {} (starting pos: {}, ending pos: {}, length: {})".format(*data))

    # here we read chromosome by chromosome to save memory and analyse them.
    N = len(gc)
    for i, chrom in enumerate(gc.chrom_names[::-1]):
        logger.info(f"==================== analysing chrom/contig {i+1}/{N} ({chrom})")
        # since we read just one contig/chromosome, the chr_list contains
        # only one contig, so we access to it with index 0

        # Performs the computation and reporting for a given chromosome
        # This call performs the analysis, and creates the HTML page
        # if HTML is creates, it fills the gc._html_list variable
        chrom_data = ChromosomeCov(gc, chrom, gc.thresholds, gc.chunksize)
        run_analysis(chrom_data, options)

        del chrom_data  # free memory for sure
        garbage.collect()

        # logging level seems to be reset to warning somewhere
        logger.setLevel(options.logging_level)

    CoverageModule(gc)

    if options.skip_multiqc is False:
        logger.info("Creating multiqc report")
        pathtocfg = sequana_data("multiqc_config.yaml", "../multiqc/")
        cmd = f"multiqc . -m sequana_coverage -f -c {pathtocfg} "

        proc = subprocess.Popen(cmd.split(), cwd=options.output_directory)
        proc.wait()

    teardown(options.output_directory)


def run_analysis(chrom, options):
    logger.info("Computing some metrics")

    if options.w_median > len(chrom) / 4:
        NW = int(len(chrom) / 4)
        if NW % 2 == 0:
            NW += 1
        logger.warning(
            "median window length is too long. \n"
            "    Setting the window length automatically to a fifth of\n"
            "    the chromosome length ({})".format(NW)
        )
    else:
        NW = options.w_median

    ######################### DEFINES OUTPUT DIR AND SAMPLE NAME  ###########
    config.output_dir = options.output_directory
    config.sample_name = os.path.basename(options.input).split(".")[0]

    directory = options.output_directory
    os.makedirs(directory, exist_ok=True)
    directory = f"{options.output_directory}/{chrom.chrom_name}"
    os.makedirs(directory, exist_ok=True)
    #########################################################################

    # compute the running median, zscore and ROIs for each chunk summarizing the
    # results in a ChromosomeCovMultiChunk instane
    logger.info("Using running median (w=%s)" % NW)
    logger.info("Number of mixture models %s " % options.k)
    results = chrom.run(
        NW, options.k, circular=options.circular, binning=options.binning, cnv_delta=options.cnv_clustering
    )
    chrom.plot_coverage(f"{directory}/coverage.png")

    if chrom.DOC < 8:
        logger.warning(
            "The depth of coverage is below 8. sequana_coverage is"
            " not optimised for such depth. You may want to "
            " increase the threshold to avoid too many false detections"
        )
    logger.info(chrom.__str__())

    # save summary and metrics
    summary = results.get_summary(caller="sequana_coverage")
    summary.to_json(f"{directory}/sequana_summary_coverage.json")

    mu = summary.data.get("fit_mu", 0)
    sigma = summary.data.get("fit_sigma", 0)
    pi = summary.data.get("pi", 0)

    logger.info(
        "Fitted central distribution (first chunk): mu=%s, sigma=%s, pi=%s"
        % (round(mu, 3), round(sigma, 3), round(pi, 3))
    )

    # some information about the ROIs found
    logger.info(
        "Searching for ROIs (threshold=[{},{}] ; double =[{},{}])".format(
            chrom.thresholds.low, chrom.thresholds.high, chrom.thresholds.low2, chrom.thresholds.high2
        )
    )
    ROIs = results.get_rois()  # results is a ChromosomeCovMultiChunk instance
    logger.info("Number of ROIs found: {}".format(len(ROIs.df)))
    logger.info("    - below average: {}".format(len(ROIs.get_low_rois())))
    logger.info("    - above average: {}".format(len(ROIs.get_high_rois())))

    # Create directory and save ROIs
    ROIs.df.to_csv(f"{directory}/rois.csv")

    logger.info(f"Creating report in {options.output_directory}. Please wait")

    if chrom._mode == "chunks":
        logger.warning("This chromosome is large. Plots in the HTML reports are skipped")

    ChromosomeCoverageModule(
        chrom,
        datatable=CoverageModule.init_roi_datatable(ROIs),
        options={"W": NW, "k": options.k, "ROIs": ROIs, "circular": options.circular},
        command=" ".join(["sequana_coverage"] + sys.argv[1:]),
    )


if __name__ == "__main__":  # pragma: no cover
    main()
