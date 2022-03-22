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
from pathlib import Path
import pandas as pd
import pysam
import numpy as np
import seaborn as sns
import shutil
import datetime
from itertools import product

from sequana import logger
from sequana.tools import reverse_complement
from sequana.gff3 import GFF3
from sequana.fasta import FastA

logger.setLevel("INFO")


class RiboDesigner(object):
    """Design probes for ribosomes depletion.

    From a complete genome assembly FASTA file and a GFF annotation file:
    - Extract genomic sequences corresponding to the selected 'seq_type'.
    - For these selected sequences, design probes of length 'probe_len' and with a space between probes of 'inter_probe_space'.
    - Detect the highest cd-hit-est identity threshold where the number of probes is inferior or equal to 'best_n_probes'.
    - Report the list of probes in bed and csv files.

    In the csv, the oligo names are in column 1 and the oligo sequences in column 2.

    :param fasta: The complete genome assembly to extract ribosome sequences from.
    :param gff: GFF annotation file of the genome assembly.
    :param output_directory: The path to the output directory.
    :param seq_type: string describing sequence annotation type (column 3 in GFF) to select rRNA from.
    :param prob_len: the size for probes in nucleotides.
    :param inter_probe_space: the space between probes in nucleotides.
    """

    def __init__(
        self,
        fasta,
        gff,
        output_directory,
        seq_type="rRNA",
        probe_len=None,
        inter_probe_space=None,
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
        self.clustered_probes_bed = self.outdir / "clustered_probes.bed"

    def get_rna_pos_from_gff(self):
        """Convert a GFF file into a pandas DataFrame filtered according to the
        self.seq_type.
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

    def get_optimal_probes(self):
        """Generate a FASTA file containing probe sequences.

        This method calculate the probe_len and inter_probe_space for each ribosomale sequence and generate FASTA probe file accordingly.

        ribo_len = probe_len * n + (inter_probe_space * (n - 1))
        <=>
        n = (ribo_len + inter_probe_space) / (prob_len + inter_probe_space)

        n being the number of probes which will be designed to have an end-to-end coverage of the ribosomale sequence.
        """
        probe_lens = range(50, 40, -1)
        inter_lens = range(20, 10, -1)

        probe_dict = {}

        with pysam.FastxFile(self.ribo_sequences_fasta) as fas:
            for seq in fas:
                seq_len = len(seq.sequence)

                for p, i in product(probe_lens, inter_lens):
                    if ((seq_len + i) / (p + i)).is_integer():
                        probe_len = p
                        inter_len = i
                        break

                for start in range(0, seq_len - probe_len + 1, probe_len + inter_len):
                    stop = start + probe_len

                    probe_dict[(seq.name, start, stop, "+")] = {
                        "probe_len": probe_len,
                        "inter_len": inter_len,
                        "sequence": f"{seq.sequence[start:stop]}",
                        "seq_id": f"{seq.name}:probe_{start}_{stop}_forward",
                    }

                    rev_seq = reverse_complement(seq.sequence)
                    add_step = int(probe_len / 2 + inter_len / 2)
                    rev_start = start + add_step
                    rev_stop = stop + add_step
                    probe_dict[(seq.name, seq_len - rev_stop, seq_len - rev_start, "-")] = {
                        "probe_len": probe_len,
                        "inter_len": inter_len,
                        "sequence": f"{rev_seq[rev_start:rev_stop]}",
                        "seq_id": f"{seq.name}:probe_{rev_start}_{rev_stop}_reverse",
                    }

        self.probes_df = pd.DataFrame(probe_dict).T
        self.probes_df.index.set_names(["chr", "start", "stop", "strand"], inplace=True)

    def export_to_fasta(self):
        """From the self.probes_df, export to FASTA and CSV files."""
        with open(self.probes_fasta, "w") as fas:
            for i, row in self.probes_df.iterrows():
                fas.write(f">{row.seq_id}\n{row.sequence}\n")

    def clustering_needed(self, force=False):

        # Do not cluster if number of probes already inferior to defined threshold
        if self.probes_df.shape[0] <= self.best_n_probes and not force:
            logger.info(f"Number of probes already inferior to {self.best_n_probes}. No clustering will be performed.")
            return False
        else:
            return True

    def cluster_probes(self, force=False):
        """Use cd-hit-est to cluster highly similar probes."""

        outdir = (
            Path(self.clustered_probes_fasta).parent
            / f"cd-hit-est-{datetime.datetime.today().isoformat(timespec='seconds', sep='_')}"
        )
        outdir.mkdir()
        log_file = outdir / "cd-hit.log"

        res_dict = {"seq_id_thres": [], "n_probes": []}
        # Add number of probes without clustering
        res_dict["seq_id_thres"].append(1)
        res_dict["n_probes"].append(self.probes_df.shape[0])

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

        if not np.isnan(best_thres):
            n_probes = df.query("seq_id_thres == @best_thres").loc[:, "n_probes"].values[0]
            logger.info(f"Best clustering threshold: {best_thres}, with {n_probes} probes.")
            shutil.copy(outdir / f"clustered_{best_thres}.fas", self.clustered_probes_fasta)
            kept_probes = [seq.name for seq in FastA(outdir / f"clustered_{best_thres}.fas")]
            self.probes_df["kept_after_clustering"] = self.probes_df.seq_id.isin(kept_probes)
            self.probes_df["bed_color"] = self.probes_df.kept_after_clustering.map(
                {True: "128,255,170", False: "255,170,128"}
            )
            self.probes_df["clustering_thres"] = best_thres
        else:
            logger.warning(f"No identity threshold was found to have as few as {self.best_n_probes} probes.")

        return df

    def export_to_csv_bed(self):
        """Export final results to CSV and BED files"""

        df = self.probes_df.query("kept_after_clustering == True")
        df.to_csv(self.clustered_probes_csv, index=False, columns=["seq_id", "sequence"])

        self.probes_df.reset_index().to_csv(
            self.clustered_probes_bed,
            index=False,
            columns=["chr", "start", "stop", "seq_id", "clustering_thres", "strand", "start", "stop", "bed_color"],
            sep="\t",
            header=False,
        )

    def fasta_to_csv(self):
        """Convert a FASTA file into a CSV file."""

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
        self.get_optimal_probes()
        self.export_to_fasta()
        if self.clustering_needed():
            self.clustered_probes_df = self.cluster_probes()
        self.export_to_csv_bed()
        # self.fasta_to_csv()
