from itertools import chain, combinations

import colorlog
import numpy as np
import pandas as pd
import rich_click as click
from jinja2 import Environment, PackageLoader
from pysam import VariantFile, VariantRecord

from sequana.gff3 import GFF3
from sequana.scripts.utils import CONTEXT_SETTINGS, common_logger

logger = colorlog.getLogger(__name__)

VCF_COLUMNS = ["chrom", "start", "stop", "type", "quality", "impact", "effect", "gene", "codon", "aa"]
NEEDED_ATB = (
    "product",
    "locus_tag",
)


def read_join_calling_vcf(vcf_file: str, qual_treshold: int):
    """Create join calling vcf dataframe"""
    with VariantFile(vcf_file) as vcf:
        # get samples information
        first = next(vcf)
        sample_names = first.samples.keys()
        variants = (
            tuple(
                [rec.chrom, rec.start, rec.stop, rec.info["TYPE"][0], rec.qual]
                + get_gene_name(rec)
                + [freq for freq in get_freqs(rec)]
            )
            for rec in chain((first,), vcf)
            if rec.qual > qual_treshold
        )
        df = pd.DataFrame(
            variants,
            columns=VCF_COLUMNS + sample_names,
        )
        return sample_names, df


def get_gene_name(rec: VariantRecord):
    """Read effect field to get information."""
    effect, info = rec.info["EFF"][0].split("(")
    impact, trans_effect, codon, aa, _, gene_name, *_ = info.split("|")
    if gene_name:
        if trans_effect:
            effect = trans_effect
        return [impact, effect, gene_name, codon, aa.split("/")[0].strip("p.")]
    return [impact, "UNKNOWN", "UNKNOWN", "NONE", "NONE"]


def get_freqs(rec):
    for info in rec.samples.values():
        try:
            yield info["AO"][0] / info["DP"]
        except TypeError:
            yield 0


@click.command(context_settings=CONTEXT_SETTINGS)
@click.option(
    "-i",
    "--input-vcf",
    type=click.Path(exists=True),
    metavar="VCF",
    nargs=1,
    required=True,
    help="Joint calling vcf using freebayes and annotated with snpEff.",
)
@click.option(
    "-g",
    "--input-gff",
    type=click.Path(exists=True),
    metavar="GFF",
    nargs=1,
    required=True,
    help="GFF used to annotate the VCF file.",
)
@click.option("-o", "--output-html", type=click.Path(), metavar="HTML", required=True, help="HTML report.")
@click.option("-t", "--title", metavar="TITLE", default="Variants comparison", help="Report title.")
@click.option(
    "-q",
    "--quality-threshold",
    type=int,
    metavar="QUAL",
    default=1,
    show_default=True,
    help="Threshold used to filter variants and keep only variants that are greater than the argument value. ",
)
@click.option(
    "-r",
    "--remove-sample",
    multiple=True,
    metavar="SAMPLE",
    help="Sample name you want to remove from cross comparisons.",
)
@click.option(
    "-s",
    "--ordered-sample",
    metavar="SAMPLE",
    help="Comma separated samples name order you want.",
)
@common_logger
def variants_comparison(
    input_vcf, input_gff, output_html, title, quality_threshold, remove_sample, ordered_sample, logger
):
    """Retrieve difference in variant content across multiple samples using joint calling vcf file."""
    sample_names, vcf = read_join_calling_vcf(input_vcf, quality_threshold)
    # remove some samples
    if remove_sample:
        vcf = vcf.drop(list(remove_sample), axis=1)
        sample_names = [s for s in sample_names if s not in remove_sample]

    if ordered_sample:
        ordered_sample = ordered_sample.split(",")
        sample_names = ordered_sample + [s for s in sample_names if s not in ordered_sample]
        vcf = vcf[VCF_COLUMNS + sample_names]

    # sometimes snpEff annotate using only locus_tag
    gff = GFF3(input_gff)
    merge_on = "gene"
    if not any(gff.df.gene.isin(vcf.gene.unique())):
        merge_on = "locus_tag"
        vcf = vcf.rename(columns={"gene": "locus_tag"})
    gff = gff.get_simplify_dataframe().filter(["gene", "locus_tag", "product"], axis=1)
    vcf = pd.merge(vcf, gff, how="left", on=merge_on)
    comparisons = {
        f"{s1}_vs_{s2}": vcf.loc[~np.isclose(vcf[s1], vcf[s2], atol=0.3)].to_dict("records")
        for s1, s2 in combinations(sample_names, r=2)
    }

    # create reporting
    template_env = Environment(loader=PackageLoader("sequana", "resources/templates/"))
    template = template_env.get_template("vcf_diff.html")
    with open(output_html, "w") as fho:
        print(
            template.render(
                title=title,
                hideable_columns=VCF_COLUMNS,
                floating_keys=set(sample_names + ["quality"]),
                comparisons=comparisons,
            ),
            file=fho,
        )
