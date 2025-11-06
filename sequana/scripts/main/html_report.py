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

import colorlog
import rich_click as click

from sequana.scripts.utils import CONTEXT_SETTINGS

logger = colorlog.getLogger(__name__)


groups = [
    {
        "name": "VCF options",
        "options": [
            "--freebayes-score",
            "--strand-ratio",
            "--frequency",
            "--forward-depth",
            "--keep-polymorphic",
            "--reverse-depth",
            "--min-depth",
            "--output-vcf-file",
            "--output-csv-file",
        ],
    },
    {
        "name": "General",
        "options": ["--output-directory"],
    },
]

click.rich_click.OPTION_GROUPS["html-report"] = groups.copy()


# not the hyphen html-report (not html_report). this is a 'bug/feature'
# @click.command(context_settings=CONTEXT_SETTINGS)
@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument("name", type=click.Path(exists=True))
@click.option("--output-directory", default=".", show_default=True)
@click.option(
    "--output-vcf-file",
    "output_vcf_file",
    required=False,
    type=click.Path(),
    show_default=True,
    default="sequana.filter.vcf",
    help="Path to output VCF file.",
)
@click.option(
    "--output-csv-file",
    "output_csv_file",
    required=False,
    type=click.Path(),
    show_default=True,
    default="sequana.filter.csv",
    help="Path to output CSV file.",
)
@click.option(
    "--freebayes-score",
    "freebayes_score",
    default=20,
    show_default=True,
    type=float,
    help="Freebayes score threshold.",
)
@click.option(
    "--strand-ratio",
    "strand_ratio",
    default=0.2,
    show_default=True,
    type=float,
    help="Minimum strand ratio.",
)
@click.option(
    "--frequency",
    "frequency",
    default=0.1,
    show_default=True,
    type=float,
    help="Minimum allele frequency.",
)
@click.option(
    "--min-depth",
    default=10,
    show_default=True,
)
@click.option(
    "--forward-depth",
    default=3,
    show_default=True,
)
@click.option(
    "--reverse-depth",
    default=3,
    show_default=True,
)
@click.option(
    "--keep-polymorphic",
    default=True,
    show_default=True,
)
def html_report(**kwargs):
    """Create a HTML report for various type of data set

    # VCF

    The VCF module takes an input VCF file and filter it before creating an HTML report.
    The options are related to the variant quality. We assume that the VCF file was
    created with sequana_variant_calling pipeline and so relied on freebayes software.
    If so, you can provide a freebayes minimal score, remove variants below a given
    frequency (--frequency), etc (please see --help for other parameters




    """
    filename = kwargs["name"]
    if filename.endswith("fastq.gz") or filename.endswith(".fastq"):
        module = "fastq"
    elif filename.endswith("vcf"):
        module = "vcf"
    else:
        logger.error("Only extensions vcf is accepted for now")
        sys.exit(1)

    if module == "vcf":
        from sequana.freebayes_vcf_filter import VCF_freebayes
        from sequana.modules_report.variant_calling import VariantCallingModule
        from sequana.utils import config

        # use the sample name based on output HTML file otherwise, user can provide one
        report_dir = kwargs["output_directory"]
        config.output_dir = report_dir
        params = {
            "freebayes_score": kwargs["freebayes_score"],
            "frequency": kwargs["frequency"],
            "min_depth": kwargs["min_depth"],
            "forward_depth": kwargs["forward_depth"],
            "reverse_depth": kwargs["reverse_depth"],
            "strand_ratio": kwargs["strand_ratio"],
            "keep_polymorphic": kwargs["keep_polymorphic"],
        }

        v = VCF_freebayes(filename)
        filter_v = v.filter_vcf(params)
        filter_v.to_csv(kwargs["output_csv_file"])
        filter_v.to_vcf(kwargs["output_vcf_file"])

        # the HTML report
        VariantCallingModule(filter_v)


# Make the script runnable
if __name__ == "__main__":
    html_report()
