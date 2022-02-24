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
import numpy as np
import seaborn as sns
import matplotlib
import shutil
import datetime

from sequana import logger
from sequana.tools import reverse_complement
from sequana.gff3 import GFF3
from sequana.fasta import FastA

logger.setLevel("INFO")


class RiboDesigner(object):
    """Documentation for RiboDesigner"""

    def __init__(
        self,
        fasta,
        gff,
        output_directory,
        seq_type="rRNA",
        probe_len=50,
        inter_probe_space=15,
        best_n_probes=384,
        force=False,
        threads=4,
    ):
        # Input
        self.fasta = fasta
        self.gff = gff
        self.seq_type = seq_type
        self.probe_len = probe_len
        self.inter_probe_space = inter_probe_space
        self.best_n_probes = best_n_probes
        self.threads = threads
        self.outdir = Path(output_directory)
        if force:
            self.outdir.mkdir(exist_ok=True)
        else:
            self.outdir.mkdir()

        # Output
        self.filtered_gff = self.outdir / "ribosome_filtered.gff"
        self.ribo_sequences_fasta = self.outdir / "ribosome_sequences.fas"
        self.probes_fasta = self.outdir / "probes_sequences.fas"
        self.clustered_probes_fasta = self.outdir / "clustered_probes.fas"
        self.clustered_probes_csv = self.outdir / "clustered_probes.csv"

    def get_rna_pos_from_gff(self):
        """Convert a GFF fil e into a pandas DataFrame filtered according to the
        seq_type.

        :param gff: GFF annotation file
        :param seq_type: string describing sequence annotation type (column 3 in GFF) to
            select rRNA from.
        """

        gff = GFF3(self.gff)
        gff.save_gff_filtered(filename=self.filtered_gff, features=[self.seq_type])
        gff_filtered = GFF3(self.filtered_gff)
        gff_filtered.to_fasta(self.fasta, self.ribo_sequences_fasta)

        genetic_types = gff.df.genetic_type.unique().tolist()
        seq_types = gff_filtered.df.ID.tolist()

        logger.info(f"Found {gff_filtered.df.shape[0]} '{self.seq_type}' entries in annotation file.")
        logger.info(f"Genetic types found in gff: {','.join(genetic_types)}")
        logger.info(f"List of '{self.seq_type}' detected: {','.join(seq_types)}")

        return gff_filtered.df

    def get_probes(self):
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
        with open(self.probes_fasta, "w") as fas_out:

            with pysam.FastxFile(self.ribo_sequences_fasta) as fas:
                for seq in fas:

                    for start in range(0, len(seq.sequence) - self.probe_len, self.probe_len + self.inter_probe_space):
                        count += 2

                        stop = start + self.probe_len

                        probe_fw = f">{seq.name}:probe_{start}_{stop}_forward\n{seq.sequence[start:stop]}\n"
                        rev_seq = reverse_complement(seq.sequence)
                        probe_rv = f">{seq.name}:probe_{start}_{stop}_reverse\n{rev_seq[start:stop]}\n"

                        fas_out.write(probe_fw + probe_rv)

        logger.info(f"{count} probes designed.")

    def cluster_probes(self):
        """Use cd-hit-est to cluster highly similar probes.

        Detect the highest cd-hit-est threshold where the number of probes is inferior or equal to best_n_probes.

        :param fasta_in: path to FASTA file with probes.
        :param fasta_out: path to FASTA file with probes after filtering.

        """

        outdir = (
            Path(self.clustered_probes_fasta).parent
            / f"cd-hit-est-{datetime.datetime.today().isoformat(timespec='seconds', sep='_')}"
        )
        outdir.mkdir()
        log_file = outdir / "cd-hit.log"

        res_dict = {"seq_id_thres": [], "n_probes": []}

        for seq_id_thres in np.arange(0.8, 1, 0.01).round(2):

            tmp_fas = outdir / f"clustered_{seq_id_thres}.fas"
            cmd = f"cd-hit-est -i {self.probes_fasta} -o {tmp_fas} -c {seq_id_thres} -n {self.threads}"
            logger.debug(f"Clustering probes with command: {cmd} (log in '{log_file}').")

            with open(log_file, "a") as f:
                subprocess.run(cmd, shell=True, check=True, stdout=f)

            res_dict["seq_id_thres"].append(seq_id_thres)
            res_dict["n_probes"].append(len(FastA(tmp_fas)))

        # Dataframe with number of probes for each cdhit identity threshold
        df = pd.DataFrame(res_dict)
        p = sns.lineplot(data=df, x="seq_id_thres", y="n_probes")
        p.axhline(self.best_n_probes, alpha=0.8, linestyle="--", color="red")

        # Extract the best identity threshold
        best_thres = df.query("n_probes <= @self.best_n_probes").seq_id_thres.max()
        print(best_thres)
        if not np.isnan(best_thres):
            n_probes = df.query("seq_id_thres == @best_thres").loc[:, "n_probes"].values[0]
            logger.info(f"Best clustering threshold: {best_thres}, with {n_probes} probes.")
            shutil.copy(outdir / f"clustered_{best_thres}.fas", self.clustered_probes_fasta)
        else:
            logger.warning(f"No identity threshold was found to have as few as {self.best_n_probes} probes.")

        return df

    def fasta_to_csv(self):
        """Convert a FASTA file into a CSV file.

        The oligo names are in column 1 and the oligo sequences in column 2.

        :param fasta_in: path to FASTA file with clustered probes.
        :param csv_out: path to CSV file with clustered probes.

        """

        seq_dict = {"id": [], "sequence": []}

        with pysam.FastxFile(self.clustered_probes_fasta) as fas:
            data = [(seq.name, seq.sequence) for seq in fas]
        logger.info(f"After clustering, keeping {len(data)} probes.")
        logger.info(f"Creating probes csv table '{self.clustered_probes_csv}'.")
        df = pd.DataFrame(data, columns=["id", "sequence"])
        df.to_csv(self.clustered_probes_csv, index=False)

        return df

    def run(self):
        self.filtered_gff_df = self.get_rna_pos_from_gff()
        self.get_probes()
        self.clustered_probes_df = self.cluster_probes()
        self.fasta_to_csv()
