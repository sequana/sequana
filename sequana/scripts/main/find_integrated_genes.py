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
import colorlog
import rich_click as click
from tqdm import tqdm
import pandas as pd
import pysam

from sequana.scripts.utils import CONTEXT_SETTINGS, common_logger
from sequana import FastQ, Cigar 

logger = colorlog.getLogger(__name__)


@click.command(context_settings=CONTEXT_SETTINGS)
@click.option(
    "--bam-file",
    help="The BAM file to introspect",
    required=True
)
@click.option(
    "--name",
    help="The name of the gene that is suppose to be integrated in the host.",
    required=True
)
@click.option(
    "--tag",
    help="""By default all output files contain the name provided with --name.
If you wish to override this behaviour you can use the --tag option. FastQ and FastA 
files will be named tag.fastq and tag.fasta""",
    required=False,
    default=None
)
@click.option(
    "--save-reads",
    help="""If provided, save reads of interest in a FastQ file named {name}.fastq and
    FastA file named {name}.fasta""",
    is_flag=True
)
@common_logger
def find_integrated_genes(**kwargs):
    """Find overlapping reads on host and possible integrated genes

    This script scans a BAM file and identifies reads that map onto a gene name provided
    by the user. Then, it checks whether these reads also map onto the other gene names.

    This is very convenient to search for an integrated plasmid or gene in a host.
    For instance a transgene in a mammal genome.

    The candidate reads are saved in various files:
    - a NAME_info.csv file with information about reads that were found of interest.
      The header is M,S,H,I,D,flag,chr,position,identifier. M, S, H, I and D indicates
      the number of bases that map (M), are soft clipped (S), hard clipped (H), inserted (I),
      or deleted (D). We also add the flag (type of mapping), the chromosome and
      position (chr, position) and finally the read identifier.
    - a NAME_reads_IDs.csv that contains the unique identifier reads
    - a FastQ if --save-reads-as-fastq option is used
    - a FastA if --save-reads-as-fasta option is used

    Note that secondary reads are ignored (flag 256, 272).

    We also removed reads that map at a single place (no overlap).

    To obtain the BAM files, you can use the sequana_mapper pipeline (https://github.com/sequana/mapper).

    The BAM file should be indexed with samtools for better performances.

    """

    gene_name = kwargs['name'].replace(" ","")
    tag = kwargs["tag"]
    if tag is None: #pragma: no cover
        tag = gene_name


    # we first read the entire dataframe
    logger.info("Scanning the BAM file. Please wait.")
    b = pysam.AlignmentFile(kwargs['bam_file'])
    results = []
    for a in tqdm(b):
        result = []
        if a.cigarstring:
            c = Cigar(a.cigarstring).as_dict()
            result.append(c.get('M', 0))
            result.append(c.get('S', 0))
            result.append(c.get('H', 0))
            result.append(c.get('I', 0))
            result.append(c.get('D', 0))
            result.append(a.flag)
            result.append(a.reference_name)
            result.append(a.pos)
            result.append(a.qname)
            results.append(result)
    df = pd.DataFrame(results)
    df.columns = ['M', 'S', 'H', 'I', 'D', 'flag', 'chr', 'position', 'identifier']

    # identify names in the BAM. could be useful in provided names are not found
    unique = df['chr'].unique()
    print(f"Found {len(unique)} unique names stored in {tag}_chr_names.txt")
    with open(f"{tag}_chr_names.txt", "w") as fout:
        for x in unique:
            fout.write(f"{x}\n")

    print(f"============== [Searching for {gene_name}] -===============")
    # identifiers reads that map on transgene
    candidates = df.query("chr == @gene_name").identifier.values
    print(f"Found {len(candidates)} hits on {gene_name}.")

    # and keep only those identifiers.
    dd = df.query("identifier in @candidates").copy()

    # ignore secondary
    dd = dd.query('flag not in [256, 272]')
    print(f"Found {len(dd)} after filtering secondary reads.")

    # another criteria is that those candidates must also map elsewhere on the host so we remove
    # identifier with only one occurence.
    dd = dd.groupby('identifier').filter(lambda x: len(x) >1)
    print(f"Found {len(dd)} reads that map in two places at least.")

    # print information on the stdout
    for chrom in dd['chr'].unique():
        #if chrom in gene_names:
        m = dd.query("chr == @chrom").position.min()
        M = dd.query("chr == @chrom").position.max()
        N = len(dd.query("chr == @chrom").position)
        if m == M: #pragma: no cover
            print(f"Found {N} hits on {chrom} chromosome at position {m}.")
        else:
            print(f"Found {N} hits on {chrom} chromosome between positions {m} and {M}.")

    print("==================[ details of found hits] =================== ")
    print(dd)

    # Saving some information.
    dd.sort_values("identifier").to_csv(f"{tag}_info.csv", sep=",", index=False)
    
    # Saving unique identifier names
    N = len(dd["identifier"].unique())
    print("==================[ Saving unique identifiers] =================== ")
    outname = f"{tag}_reads_IDs.csv"
    print(f"Number of unique identifiers: {N}. Stored in {outname}")
    pd.Series(dd["identifier"].unique()).to_csv(outname, sep=",", index=False, header=None)


    if gene_name not in dd['chr'].unique(): #pragma: no cover
        print(f"Not hits were found for {gene_name}. Look at the chr_name.txt file to check your gene's name")

    # rescan and save FastQ file
    identifiers = list(dd['identifier'].values)
    if kwargs['save_reads']:
        print("==================[ Saving reads] =================== ")
        outfastq = f"{tag}.fastq"
        outfasta = f"{tag}.fasta"

        logger.info(f"Saving reads of interest in {outfastq} and {outfasta}")
        b = pysam.AlignmentFile(kwargs['bam_file'])
        with open(outfastq, "w") as fastq, open(outfasta, "w") as fasta:
            for a in tqdm(b):
                if a.qname and a.qname in identifiers:
                    fastq.write(f"@{a.qname}\n{a.query}\n+\n{a.qual}\n")
                    fasta.write(f">{a.qname}\n{a.query}\n")
                    # we pop up identifiers to avoid writting reads several times
                    identifiers.remove(a.qname)
                if len(identifiers) == 0:
                    break






