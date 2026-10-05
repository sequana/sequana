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

from sequana.fasta import FastAComparison
from sequana.scripts.utils import CONTEXT_SETTINGS, common_logger

logger = colorlog.getLogger(__name__)


@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument("filename1", type=click.Path(exists=True))
@click.argument("filename2", type=click.Path(exists=True))
@click.option(
    "--rc-aware",
    is_flag=True,
    help="consider a sequence and its reverse complement as identical",
)
@click.option(
    "--strict-case",
    is_flag=True,
    help="case-sensitive comparison. A soft-masked sequence then differs from its unmasked version",
)
@click.option("--ignore-gaps", is_flag=True, help="drop the '-' and '.' alignment characters before comparing")
@click.option("-o", "--output", help="save the comparison as a TSV file")
@common_logger
def fasta_compare(**kwargs):
    """Compare the sequences of two FastA files independently of their names

    Sequences are paired using the MD5 checksum of their content, so that
    identical sequences are matched even when their identifiers differ (e.g. *1*
    in one file versus *NC_045512.2* in the other one). Unmatched sequences are
    reported as orphans on each side.

    \b
        sequana fasta-compare reference.fa assembly.fa
        sequana fasta-compare reference.fa assembly.fa --rc-aware -o report.tsv

    The exit code is 0 when both files have the same sequence content and 1
    otherwise, so that the command can be used in a script.

    """
    c = FastAComparison(
        kwargs["filename1"],
        kwargs["filename2"],
        ignore_case=not kwargs["strict_case"],
        ignore_gaps=kwargs["ignore_gaps"],
        rc_aware=kwargs["rc_aware"],
    )

    matches = c.matches
    orphans1 = c.orphans1
    orphans2 = c.orphans2
    N1 = sum(len(x) for x in c.data1.values())
    N2 = sum(len(x) for x in c.data2.values())
    Nmatched1 = N1 - len(orphans1)
    Nmatched2 = N2 - len(orphans2)

    print(f"# {kwargs['filename1']}: {N1} sequences, {Nmatched1} matched")
    print(f"# {kwargs['filename2']}: {N2} sequences, {Nmatched2} matched")

    if matches:
        print("\n# matches")
        for names1, names2, length, rc in matches:
            tag = " (reverse complement)" if rc else ""
            print(f"{','.join(names1)} == {','.join(names2)}\t{length} bp{tag}")

    for filename, orphans in [(kwargs["filename1"], orphans1), (kwargs["filename2"], orphans2)]:
        if orphans:
            print(f"\n# orphans in {filename}")
            for name, length in orphans:
                print(f"{name}\t{length} bp")

    for filename, duplicates in [(kwargs["filename1"], c.duplicates1), (kwargs["filename2"], c.duplicates2)]:
        if duplicates:
            print(f"\n# duplicated sequences in {filename}")
            for names in duplicates:
                print(",".join(names))

    if kwargs["output"]:
        with open(kwargs["output"], "w") as fout:
            fout.write("status\tnames1\tnames2\tlength\treverse_complement\n")
            for names1, names2, length, rc in matches:
                fout.write(f"match\t{','.join(names1)}\t{','.join(names2)}\t{length}\t{rc}\n")
            for name, length in orphans1:
                fout.write(f"orphan1\t{name}\t\t{length}\t\n")
            for name, length in orphans2:
                fout.write(f"orphan2\t\t{name}\t{length}\t\n")
        logger.info(f"Saved comparison in {kwargs['output']}")

    if c.identical:
        logger.info("The two files contain the same sequences.")
    else:
        logger.warning("The two files differ.")
        sys.exit(1)
