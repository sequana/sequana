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
import sys
import json
import subprocess
from pathlib import Path

import pandas as pd
import pylab
import pysam
import numpy as np
import seaborn as sns
import shutil
import datetime
from itertools import product

from sequana import logger
from sequana.tools import reverse_complement
from sequana.fasta import FastA

logger.setLevel("INFO")


class RiboDesigner(object):
    """Design probes for ribosomes depletion.

    From a complete genome assembly FASTA file and a GFF annotation file:

    - Extract genomic sequences corresponding to the selected ``seq_type``.
    - For these selected sequences, design probes computing probe length and inter probe space according to the length of the ribosomale sequence.
    - Detect the highest cd-hit-est identity threshold where the number of probes is inferior or equal to ``max_n_probes``.
    - Report the list of probes in BED and CSV files.

    In the CSV, the oligo names are in column 1 and the oligo sequences in column 2.

    :param fasta: The FASTA file with complete genome assembly to extract ribosome sequences from.
    :param gff: GFF annotation file of the genome assembly.
    :param output_directory: The path to the output directory.
    :param seq_type: string describing sequence annotation type (column 3 in GFF) to select rRNA from.
    :param max_n_probes: Max number of probes to design
    :param force:  If the `output_directory` already exists, overwrite it.
    :param threads: Number of threads to use in cd-hit clustering.
    :param float identity_step: step to scan the sequence identity (between 0 and 1) defaults to 0.01.
    """

    def __init__(
        self,
        fasta,
        gff,
        output_directory,
        seq_type="rRNA",
        max_n_probes=384,
        force=False,
        threads=4,
        identity_step=0.01,
        **kwargs
    ):
        # Input
        self.fasta = fasta
        self.gff = gff
        self.seq_type = seq_type
        self.max_n_probes = max_n_probes
        self.threads = threads
        self.outdir = Path(output_directory)
        self.identity_step = identity_step

        if force:
            self.outdir.mkdir(exist_ok=True)
        else:
            try:
                self.outdir.mkdir()
            except FileExistsError as err:
                logger.error(f"Output directory {output_directory} exists. Use --force or set force=True")
                sys.exit(1)

        # Output
        self.filtered_gff = self.outdir / "ribosome_filtered.gff"
        self.ribo_sequences_fasta = self.outdir / "ribosome_sequences.fas"
        self.probes_fasta = self.outdir / "probes_sequences.fas"
        self.clustered_probes_fasta = self.outdir / "clustered_probes.fas"
        self.clustered_probes_csv = self.outdir / "clustered_probes.csv"
        self.clustered_probes_bed = self.outdir / "clustered_probes.bed"
        self.output_json = self.outdir / "ribodesigner.json"

        self.json = {"max_n_probes": max_n_probes, "identity_step": identity_step, "feature": seq_type}

    def get_rna_pos_from_gff(self):
        """Convert a GFF file into a pandas DataFrame filtered according to the
        self.seq_type.
        """

        gff = pd.read_csv(
            self.gff,
            sep="\t",
            comment="#",
            names=["seqid", "source", "seq_type", "start", "end", "score", "strand", "phase", "attributes"],
        )

        filtered_gff = gff.query("seq_type == @self.seq_type")

        with pysam.Fastafile(self.fasta) as fas:
            with open(self.ribo_sequences_fasta, "w") as fas_out:

                for row in filtered_gff.itertuples():
                    region = f"{row.seqid}:{row.start}-{row.end}"
                    seq_record = f">{region}\n{fas.fetch(region=region)}\n"
                    fas_out.write(seq_record)

        seq_types = gff.seq_type.unique().tolist()

        logger.info(f"Genetic types found in gff: {','.join(seq_types)}")
        logger.info(f"Found {filtered_gff.shape[0]} '{self.seq_type}' entries in the annotation file.")
        logger.debug(f"\t" + filtered_gff.to_string().replace("\n", "\n\t"))

        filtered_gff.to_csv(self.filtered_gff)
        
    def _get_probe_and_step_len(self, seq):
        """Calculates the probe_len and inter_probe_space for a ribosomal sequence.

        ribo_len = probe_len * n + (inter_probe_space * (n - 1))
        <=>
        n = (ribo_len + inter_probe_space) / (prob_len + inter_probe_space)
        """
        seq_len = len(seq.sequence)

        probe_lens = range(50, 40, -1)
        inter_probe_space = range(20, 10, -1)

        for probe_len, inter_probe_space in product(probe_lens, inter_probe_space):
            if ((seq_len + inter_probe_space) / (probe_len + inter_probe_space)).is_integer():
                return probe_len, inter_probe_space

        raise ValueError(f"No correct probe length/inter probe space combination was found for {seq.name}")

    def _get_probes_df(self, seq, probe_len, step_len):
        """Generate the Dataframe with probes information.

        Design probes to have end-to-end coverage on the + strand and fill the inter_probe_space present on the + strand with probes designed on the - strand.

        :param seq: A pysam sequence object.
        :param prob_len: The length of the probes calculated by self._get_probe_and_step_len.
        :param step_len: The length of the inter-probe space calculated by self._get_probe_and_step_len.
        :param strand: The strand on which probes are designed.
        """

        # + strand probes
        starts = [start for start in range(0, len(seq.sequence) - probe_len + 1, probe_len + step_len)]
        stops = [start + probe_len for start in starts]

        df = pd.DataFrame({"name": seq.name, "start": starts, "stop": stops, "strand": "+", "score": 0})
        df["sequence"] = [seq.sequence[row.start : row.stop] for row in df.itertuples()]
        df["seq_id"] = df["name"] + f"_+_" + df["start"].astype(str) + "_" + df["stop"].astype(str)

        # - strand probes
        sequence = reverse_complement(seq.sequence)
        # Starts reverse probes to be centered on inter_probe_space of the forward probes
        rev_starts = [int((starts[i + 1] + starts[i]) / 2) for i in range(0, len(starts) - 1)]
        rev_stops = [start + probe_len for start in rev_starts]

        df_rev = pd.DataFrame({"name": seq.name, "start": rev_starts, "stop": rev_stops, "strand": "-", "score": 0})
        df_rev["sequence"] = [sequence[row.start : row.stop] for row in df_rev.itertuples()]
        df_rev["seq_id"] = df_rev["name"] + f"_-_" + df_rev["start"].astype(str) + "_" + df_rev["stop"].astype(str)

        # Transform to bed coordinates for the reverse_complement
        df_rev["start"] = len(sequence) - df_rev["start"]
        df_rev["stop"] = len(sequence) - df_rev["stop"]
        df_rev.rename(columns={"start": "stop", "stop": "start"}, inplace=True)

        return pd.concat([df, df_rev])

    def get_all_probes(self):
        """Run all probe design and concatenate results in a single DataFrame."""

        probes_dfs = []

        with pysam.FastxFile(self.ribo_sequences_fasta) as fas:
            for seq in fas:

                probe_len, step_len = self._get_probe_and_step_len(seq)
                probes_dfs.append(self._get_probes_df(seq, probe_len, step_len))

        self.probes_df = pd.concat(probes_dfs)
        self.probes_df["kept_after_clustering"] = True
        self.probes_df["bed_color"] = self.probes_df.kept_after_clustering.map(
            {True: "21,128,0", False: "128,64,0"}
        )

    def export_to_fasta(self):
        """From the self.probes_df, export to FASTA and CSV files."""

        with open(self.probes_fasta, "w") as fas:
            for row in self.probes_df.itertuples():
                fas.write(f">{row.seq_id}\n{row.sequence}\n")

    def clustering_needed(self, force=False):
        """Checks if a clustering is needed.

        :param force: force clustering even if unecessary.
        """

        # Do not cluster if number of probes already inferior to defined threshold
        if not force and self.probes_df.shape[0] <= self.max_n_probes:
            logger.info(f"Number of probes already inferior to {self.max_n_probes}. No clustering will be performed.")
            return False
        else:
            return True

    def cluster_probes(self):
        """Use cd-hit-est to cluster highly similar probes."""

        outdir = (
            Path(self.clustered_probes_fasta).parent
            / f"cd-hit-est-{datetime.datetime.today().isoformat(timespec='seconds', sep='_')}"
        )
        outdir.mkdir()
        log_file = outdir / "cd-hit.log"

        res_dict = {"seq_id_thres": [], "n_probes": []}

        for seq_id_thres in np.arange(0.8, 1, self.identity_step).round(3):

            tmp_fas = outdir / f"clustered_{seq_id_thres}.fas"
            cmd = f"cd-hit-est -i {self.probes_fasta} -o {tmp_fas} -c {seq_id_thres} -n {self.threads}"
            logger.debug(f"Clustering probes with command: {cmd} (log in '{log_file}').")

            with open(log_file, "a") as f:
                subprocess.run(cmd, shell=True, check=True, stdout=f)

            res_dict["seq_id_thres"].append(seq_id_thres)
            res_dict["n_probes"].append(len(FastA(tmp_fas)))

        # Add number of probes without clustering
        res_dict["seq_id_thres"].append(1)
        res_dict["n_probes"].append(self.probes_df.shape[0])

        self.json["results"] = res_dict

        # Dataframe with number of probes for each cdhit identity threshold
        pylab.clf()
        df = pd.DataFrame(res_dict)
        p = sns.lineplot(data=df, x="seq_id_thres", y="n_probes", markers=["o"])
        p.axhline(self.max_n_probes, alpha=0.8, linestyle="--", color="red", label="max number of probes requested")
        pylab.xlabel("Sequence identity", fontsize=16)
        pylab.ylabel("Number of probes", fontsize=16)


        # Extract the best identity threshold
        best_thres = df.query("n_probes <= @self.max_n_probes").seq_id_thres.max()

        if not np.isnan(best_thres):
            n_probes = df.query("seq_id_thres == @best_thres").loc[:, "n_probes"].values[0]

            self.json['n_probes'] = int(n_probes)
            self.json['best_thres'] = best_thres

            logger.info(f"Best clustering threshold: {best_thres}, with {n_probes} probes.")
            shutil.copy(outdir / f"clustered_{best_thres}.fas", self.clustered_probes_fasta)
            kept_probes = [seq.name for seq in FastA(outdir / f"clustered_{best_thres}.fas")]
            self.probes_df["kept_after_clustering"] = self.probes_df.seq_id.isin(kept_probes)
            self.probes_df["bed_color"] = self.probes_df.kept_after_clustering.map(
                {True: "21,128,0", False: "128,64,0"}
            )
            self.probes_df["clustering_thres"] = best_thres
            pylab.plot(best_thres, n_probes, 'o', label="Final number of probes")
            pylab.legend()
        else:
            logger.warning(f"No identity threshold was found to have as few as {self.max_n_probes} probes.")

        self.clustering_df = df.sort_values("seq_id_thres")
        
        return self.probes_df.query("kept_after_clustering == True")

    def export_to_csv_bed(self):
        """Export final results to CSV and BED files"""

        if self.clustering_needed():
            df = self.cluster_probes()
        else:
            df = self.probes_df
            
        df.to_csv(self.clustered_probes_csv, index=False, columns=["seq_id", "sequence"])

        self.probes_df.to_csv(
            self.clustered_probes_bed,
            sep="\t",
            index=False,
            header=None,
            columns=["name", "start", "stop", "sequence", "score", "strand", "start", "stop", "bed_color"],
        )

    def export_to_json(self):
        with open(self.output_json, "w") as fout:
            json.dump(self.json, fout, indent=4, sort_keys=True)

    def run(self):
        self.get_rna_pos_from_gff()
        self.get_all_probes()
        self.export_to_fasta()
        self.export_to_csv_bed()
        self.export_to_json()
