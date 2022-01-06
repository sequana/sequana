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
"""Ribodesigner module"""

import subprocess
import sys
from pathlib import Path
import tempfile as tmp

import pandas as pd
import pysam

from sequana import logger
from sequana.tools import reverse_complement
from sequana.gff3 import GFF3

logger.setLevel("INFO")


def get_rna_pos_from_gff(gff, gff_filtered, ref_fasta, fasta_out, seq_type="rRNA"):
    """Convert a GFF file into a pandas DataFrame filtered according to the
    seq_type.

    :param gff: GFF annotation file
    :param seq_type: string describing sequence annotation type (column 3 in GFF) to
        select rRNA from.
    """

    gff = GFF3(gff)
    gff.save_gff_filtered(filename=gff_filtered, features=[seq_type])
    gff_filtered = GFF3(gff_filtered)
    gff_filtered.to_fasta(ref_fasta, fasta_out)

    logger.info(f"Found {gff_filtered.df.shape[0]} '{seq_type}' entries in annotation file.")

    return gff_filtered.df


def get_probes(fasta_in, fasta_out, probe_len=50, inter_probe_space=15):
    """Generate a FASTA file containing probe sequences.

    From a FASTA file of nucleotide sequence to be targetted, generate a FASTA
    file with probe designed of size 'probe_len' and with a space between probes
    of 'inter_probe_space'

    :param fasta_in: path to FASTA file with extracted rRNA regions.
    :param fasta_out: path to FASTA file where probes will be exported to.
    :param prob_len: the size for probes in nucleotides.
    :param inter_probe_space: the space between probes in nucleotides.
    """

    count = 0
    with open(fasta_out, "w") as fas_out:

        with pysam.FastxFile(fasta_in) as fas:
            for seq in fas:

                for start in range(0, len(seq.sequence) - probe_len, probe_len + inter_probe_space):
                    count += 2

                    stop = start + probe_len

                    probe_fw = f">{seq.name}:probe_{start}_{stop}_forward\n{seq.sequence[start:stop]}\n"
                    rev_seq = reverse_complement(seq.sequence)
                    probe_rv = f">{seq.name}:probe_{start}_{stop}_reverse\n{rev_seq[start:stop]}\n"

                    fas_out.write(probe_fw + probe_rv)

    logger.info(f"{count} probes designed.")


def cluster_probes(fasta_in, fasta_out, seq_id_thres=0.80, threads=4):
    """Use cd-hit-est to cluster highly similar probes

    :param fasta_in: path to FASTA file with probes.
    :param fasta_out: path to FASTA file with probes after filtering.
    """
    log_file = Path(fasta_out).with_suffix(".log")

    cmd = f"cd-hit-est -i {fasta_in} -o {fasta_out} -c {seq_id_thres} -n {threads}"
    logger.info(f"Clustering probes with command: {cmd} (log in '{log_file}').")

    with open(log_file, "w") as f:
        subprocess.run(cmd, shell=True, check=True, stdout=f)


def fasta_to_csv(fasta_in, csv_out):
    """Convert a FASTA file into a CSV file.

    The oligo names are in column 1 and the oligo sequences in column 2.

    :param fasta_in: path to FASTA file with clustered probes.
    :param csv_out: path to CSV file with clustered probes.

    """

    seq_dict = {"id": [], "sequence": []}

    with pysam.FastxFile(fasta_in) as fas:
        data = [(seq.name, seq.sequence) for seq in fas]
    logger.info(f"After clustering, keeping {len(data)} probes.")
    logger.info(f"Creating probes csv table '{csv_out}'.")
    pd.DataFrame(data, columns=["id", "sequence"]).to_csv(csv_out, index=False)
