# sequana/cli/fastq_split.py
import rich_click as click

from sequana import logger
from sequana.fastq_splitter import FastqSplitter


@click.command(name="fastq-split")
@click.argument("input_fastq", type=click.Path(exists=True))
@click.option("--by-size", type=int, help="Split file into chunks of N reads each.")
@click.option("--by-part", type=int, help="Split file into N equal parts.")
@click.option(
    "--pattern",
    default="output.#.fastq.gz",
    show_default=True,
    help="Output filename pattern (use # as placeholder for chunk number).",
)
@click.option("--gzip", "gzip_output", is_flag=True, help="Compress output with gzip.")
@click.option(
    "--buffer-size",
    type=int,
    default=10000,
    show_default=True,
    help="Number of reads to buffer before writing to disk.",
)
def fastq_split(input_fastq, by_size, by_part, pattern, gzip_output, buffer_size):
    """
    Split a FASTQ file into smaller parts (by number of reads or by number of parts).

    Examples:

        sequana fastq-split reads.fastq --by-size 1000000 --pattern chunk.#.fastq.gz

        sequana fastq-split reads.fastq.gz --by-part 10 --pattern sample.#.fastq.gz
    """
    if (by_size and by_part) or not (by_size or by_part):
        raise click.UsageError("Use either --by-size or --by-part, not both.")

    splitter = FastqSplitter(input_fastq)

    if by_size:
        logger.info(f"Splitting {input_fastq} into chunks of {by_size} reads each.")
        splitter.by_size(by_size, pattern, gzip_output, buffer_size)
    else:
        logger.info(f"Splitting {input_fastq} into {by_part} equal parts.")
        splitter.by_part(by_part, pattern, gzip_output, buffer_size)

    logger.info("Splitting completed successfully.")
