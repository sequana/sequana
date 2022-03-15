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

import click
import colorlog

from .utils import CONTEXT_SETTINGS


logger = colorlog.getLogger(__name__)


@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument("name", type=click.Path(exists=True), nargs=-1)
@click.option(
    "--module",
    required=False,
    type=click.Choice(["bamqc", "bam", "fasta", "fastq", "gff", "vcf"]),
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

            BAMQCModule(name, "bamqc.html")
    elif module == "fasta":  # there is no module per se. HEre we just call FastA.summary()
        from sequana.fasta import FastA

        for name in names:
            f = FastA(name)
            f.summary()
    elif module == "fastq":  # there is no module per se. HEre we just call FastA.summary()
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
            print(f"#filename: {filename}")
            print("#Number of entries per genetic type:")
            print(ff.df.value_counts("genetic_type").to_string())
            print("#Number of duplicated attribute (if any) per attribute:")
            ff.get_duplicated_attributes_per_genetic_type()
    elif module == "vcf":
        from sequana.freebayes_vcf_filter import VCF_freebayes

        for filename in names:
            print(f"#filename: {filename}")
            vcf = VCF_freebayes(filename)
            columns = (
                "chr",
                "position",
                "depth",
                "reference",
                "alternative",
                "freebayes_score",
                "strand_balance",
                "frequency",
            )
            print(",".join(columns))
            for variant in vcf.get_variants():
                resume = variant.resume
                print(",".join([str(resume[col]) for col in columns]))