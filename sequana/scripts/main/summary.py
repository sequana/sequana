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
import json
import sys
from datetime import datetime
from pathlib import Path

import colorlog
import rich_click as click
from rich.console import Console
from rich.table import Table

from sequana.scripts.utils import CONTEXT_SETTINGS

logger = colorlog.getLogger(__name__)
console = Console()


def _convert_to_serializable(obj):
    """Convert numpy/pandas types to native Python types for JSON."""
    import numpy as np

    if isinstance(obj, dict):
        return {k: _convert_to_serializable(v) for k, v in obj.items()}
    elif isinstance(obj, (list, tuple)):
        return [_convert_to_serializable(v) for v in obj]
    elif isinstance(obj, (np.integer, np.floating)):
        return obj.item()
    return obj


def _make_json_output(filename, format_type, stats):
    result = {
        "filename": str(filename),
        "format": format_type,
        "timestamp": datetime.utcnow().isoformat() + "Z",
        "stats": _convert_to_serializable(stats),
    }
    return result


def _print_table(filename, format_type, stats):
    """Print stats in formatted rich table."""
    table = Table(title=f"{format_type.upper()} Summary: {Path(filename).name}")
    table.add_column("Metric", style="cyan")
    table.add_column("Value", style="green")

    for key, value in stats.items():
        if isinstance(value, dict):
            table.add_row(key, "")
            for k, v in value.items():
                table.add_row(f"  {k}", str(v))
        else:
            if isinstance(value, float):
                table.add_row(key, f"{value:.2f}")
            else:
                table.add_row(key, str(value))

    console.print(table)


def _print_fasta_contigs_table(contigs):
    """Print per-contig info for FASTA (limit to 200)."""
    if not contigs:
        return

    table = Table(title="Top Contigs by Length (max 200)")
    table.add_column("Contig Name", style="cyan")
    table.add_column("Length", style="green")
    table.add_column("GC%", style="yellow")
    table.add_column("A", style="magenta")
    table.add_column("C", style="magenta")
    table.add_column("G", style="magenta")
    table.add_column("T", style="magenta")
    table.add_column("N", style="magenta")

    for contig_info in contigs[:200]:
        table.add_row(
            contig_info["name"],
            str(contig_info["length"]),
            f"{contig_info['gc']:.1f}",
            str(contig_info["A"]),
            str(contig_info["C"]),
            str(contig_info["G"]),
            str(contig_info["T"]),
            str(contig_info["N"]),
        )

    console.print(table)


@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument("name", type=click.Path(exists=True), nargs=-1)
@click.option(
    "--module",
    required=False,
    type=click.Choice(["bamqc", "bam", "fna", "fasta", "fastq", "gff", "vcf", "sam"]),
)
@click.option("--output-file", required=False, type=click.Path())
@click.option("--output-json", required=False, type=click.Path(), help="Export stats to JSON file")
def summary(**kwargs):
    """Create a HTML report for various types of NGS formats.

    Supported modules:

    - bamqc
    - fastq

    This processes all files in the given pattern (in back-quotes)
    sequentially and produces one HTML file per input.

    Other modules work the same way. For example, for FastQ files::

        sequana summary one_input.fastq
        sequana summary `ls *fastq`

    Export to JSON::

        sequana summary input.fastq --output-json stats.json
    """
    names = kwargs["name"]
    module = kwargs["module"]
    output_json = kwargs["output_json"]

    results = []

    if module is None:
        if names[0].endswith("fastq.gz") or names[0].endswith(".fastq"):
            module = "fastq"
        elif names[0].endswith(".bam"):
            module = "bam"
        elif names[0].endswith(".sam"):
            module = "bam"
        elif names[0].endswith(".gff") or names[0].endswith("gff3"):
            module = "gff"
        elif names[0].endswith("fasta.gz") or names[0].endswith(".fasta"):
            module = "fasta"
        elif names[0].endswith("fa.gz") or names[0].endswith(".fa"):
            module = "fasta"
        elif names[0].endswith("fna.gz") or names[0].endswith(".fna"):
            module = "fasta"
        elif names[0].endswith(".vcf") or names[0].endswith(".vcf.gz") or names[0].endswith(".bcf"):
            module = "vcf"
        else:
            logger.error(
                "Only extensions fastq, fasta, bam, sam, gff, gff3, vcf, vcf.gz and bcf are recognised. please use --module to tell us about the type of the input files"
            )
            sys.exit(1)

    if module == "bamqc":
        for name in names:
            print(f"Processing {name}")
            from sequana.modules_report.bamqc import BAMQCModule

            BAMQCModule(name, "bamqc.html")
    elif module == "fasta":
        from sequana.fasta import FastA

        for name in names:
            f = FastA(name)
            stats = f.get_stats()

            # Calculate GC content for overall stats
            gc_overall = f.GC_content()
            stats["GC_content"] = gc_overall

            # Collect per-contig info sorted by length
            contigs = []
            for i, contig_name in enumerate(f.names):
                seq = f._fasta.fetch(contig_name)
                length = f.lengths[i]
                a_count = seq.count("A") + seq.count("a")
                c_count = seq.count("C") + seq.count("c")
                g_count = seq.count("G") + seq.count("g")
                t_count = seq.count("T") + seq.count("t")
                n_count = seq.count("N") + seq.count("n")
                gc_denom = a_count + c_count + g_count + t_count
                gc = (g_count + c_count) / gc_denom * 100 if gc_denom else 0

                contigs.append(
                    {
                        "name": contig_name,
                        "length": length,
                        "A": a_count,
                        "C": c_count,
                        "G": g_count,
                        "T": t_count,
                        "N": n_count,
                        "gc": gc,
                    }
                )

            # Sort by length descending
            contigs.sort(key=lambda x: x["length"], reverse=True)

            if output_json:
                results.append(_make_json_output(name, "fasta", stats))
            else:
                _print_table(name, "fasta", stats)
                _print_fasta_contigs_table(contigs)
    elif module == "fastq":
        from sequana import FastQC

        for filename in names:
            ff = FastQC(filename, max_sample=1e6, verbose=False)
            stats_df = ff.get_stats()
            if output_json:
                stats_dict = stats_df.iloc[0].to_dict()
                results.append(_make_json_output(filename, "fastq", stats_dict))
            else:
                stats_dict = stats_df.iloc[0].to_dict()
                _print_table(filename, "fastq", stats_dict)
    elif module == "bam":
        from sequana import BAM

        for filename in names:
            ff = BAM(filename)
            stats = ff.get_stats()
            if output_json:
                results.append(_make_json_output(filename, "bam", stats))
            else:
                _print_table(filename, "bam", stats)
    elif module == "sam":
        from sequana import SAM

        for filename in names:
            ff = SAM(filename)
            stats = ff.get_stats()
            if output_json:
                results.append(_make_json_output(filename, "sam", stats))
            else:
                _print_table(filename, "sam", stats)
    elif module == "gff":
        from sequana import GFF3

        for filename in names:
            ff = GFF3(filename)
            feature_counts = ff.df["genetic_type"].value_counts().to_dict()
            stats = {"feature_counts": feature_counts}
            if output_json:
                results.append(_make_json_output(filename, "gff", stats))
            else:
                _print_table(filename, "gff", stats)
    elif module == "vcf":
        from sequana.variants import VariantFile

        for filename in names:
            vcf = VariantFile(filename, progress=True)
            df = vcf.df
            variant_counts = df["type"].value_counts().to_dict()
            stats = {"variant_counts": variant_counts}
            if output_json:
                results.append(_make_json_output(filename, "vcf", stats))
            else:
                _print_table(filename, "vcf", stats)
                if kwargs["output_file"]:
                    df.to_csv(kwargs["output_file"], index=False)

    if output_json and results:
        with open(output_json, "w") as fh:
            json.dump(results, fh, indent=2)
