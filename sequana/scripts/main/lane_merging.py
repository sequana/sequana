#  This file is part of Sequana software
#
#  Copyright (c) 2016-2022 - Sequana Development Team
#
#  Distributed under the terms of the 3-clause BSD license.
#  The full license is in the LICENSE file, distributed with this software.
#
#  website: https://github.com/sequana/sequana
#  documentation: http://sequana.readthedocs.io
#
##############################################################################
import os
import glob
import re
import subprocess

from shutil import which
from subprocess import STDOUT

from snakemake import shell as shellcmd

import click
import colorlog

from .utils import CONTEXT_SETTINGS, OptionEatAll


logger = colorlog.getLogger(__name__)


class Common(object):
    def __init__(self, pattern, lanes):
        """

        :param pattern: example */*fastq.gz
        :param lanes: lanes to keep (and merge)

        """
        # sanity checks.
        assert len(lanes) >= 2, "to merge lanes, you need at least 2!"
        assert min(lanes) > 0, "lanes must be positive numbers"
        assert len(set(lanes)) == len(lanes), "lanes number must be unique"
        assert max(lanes) <= 8, "lanes must be < 8 "

        self.pattern = pattern
        self.lanes = ["L00{}".format(lane) for lane in lanes]
        self.Nlanes = len(lanes)
        self.discover_files()

    def discover_files(self):
        """Extract unique sample names"""
        self.filenames = glob.glob(os.path.abspath(self.pattern))
        print("Found {} files ".format(len(self.filenames)))

        if len(self.filenames) == 0:
            ValueError("No files found. Are you in the correct path ? ")

        print("Identifying lanes from input filenames...")
        lanes = str(self.lanes).replace(" ", "")
        self.filenames = [x for x in self.filenames if re.search("_L00{}_".format(lanes), x)]

        if len(self.filenames) == 0:
            ValueError("No files found with L00 tag. Are you in the correct path ? ")

        self.sampleIDs = sorted(list(set([x.split("_L00")[0] for x in self.filenames])))

        assert len(self.sampleIDs) == len(set(self.sampleIDs))


class LaneMerger(Common):
    """
    In backspace, the sample sheet contains the project, samplesID and name.

    The files are stored in a tree structure as follows::

        <Projet>/<sampleID>/<name>_R1_.fastq.gz
        <Projet>/<sampleID>/<name>_R2_.fastq.gz

    This script works at the sampleID level. Even though there are found
    in subsirectories, at the end we store everything in a flat structure.

    """

    def __init__(self, pattern="*/*fastq.gz", outdir="merging", threads=4, queue=None, lanes=[], force=False):
        super(LaneMerger, self).__init__(pattern, lanes)

        self._outdir = None
        self.outdir = outdir
        self.threads = threads
        self.queue = queue
        self.force = force

        print("Number of samples: {}".format(len(self.sampleIDs)))

    def _set_outdir(self, outdir):
        self._outdir = outdir
        self._create_outdir()

    def _get_outdir(self):
        return self._outdir

    outdir = property(_get_outdir, _set_outdir)

    def _create_outdir(self):
        if os.path.exists(self.outdir):
            pass
        else:
            print("Creating ./{} directory".format(self.outdir))
            os.makedirs(self.outdir)

    def is_paired(self, sampleID):
        filenames = [x for x in self.filenames if x.startswith(sampleID + "_L00")]
        R1 = sum([1 for filename in filenames if "_R1_" in filename])
        R2 = sum([1 for filename in filenames if "_R2_" in filename])

        if R1 == self.Nlanes and R2 == self.Nlanes:
            return True
        elif R1 == self.Nlanes and R2 == 0:
            return False
        else:
            msg = "Sample {} issue. Found {} R1 and {} R2. "
            msg += "Expected {} and {} for paired data and {} and 0 for single read"
            msg = msg.format(name, R1, R2, self.Nlanes, self.Nlanes, self.Nlanes)
            raise ValueError(msg)

    def get_pigz_cmd(self, sampleID, RX):
        assert RX in ["R1", "R2"]
        params = {"thread": self.threads, "sampleID": sampleID, "RX": RX}

        names = [x for x in self.filenames if x.startswith(sampleID) and RX in x]

        print("  Found {} files for sample {} ({})".format(len(names), sampleID, RX))
        msg = "For sample ID {}, found non unique sample names {} ?".format(sampleID, names)
        assert len(set(names)) == self.Nlanes, msg
        params["name"] = sampleID.split("/")[-1]
        params["outdir"] = self.outdir

        # IMPORTANT TO SORT the filenames so that L1 follows L2 for R1 and R2
        # files
        params["filenames"] = " ".join(sorted(names))

        output = "{outdir}/{name}_{RX}_001.fastq".format(**params)

        if os.path.exists(output + ".gz") and self.force is False:
            raise IOError("{} exists already".format(output))
        # Here, we should get 4 files with identical NAME (prefix), which should
        # be used for the output.
        # 53_S53_L001_R2_001.fastq.gz

        cmd = "#!/bin/sh"
        cmd += "\nunpigz -p {thread} -c {filenames} > " + output
        cmd += "\npigz --force -p {thread} " + output
        cmd = cmd.format(**params)
        return cmd

    def run(self, dry_run=False):

        processes = []
        for name in self.sampleIDs:
            print("{}, paired is {}".format(name, self.is_paired(name)))

            # R1. Note the usage of wrap using --wrap " your command"
            if self.queue == "biomicspole":
                sbatch_command = "sbatch -c {thread} -A biomics --qos biomics -p biomics"
            elif self.queue == "common":
                sbatch_command = "sbatch -c {thread} --qos fast"
            else:
                sbatch_command = "sbatch -c {thread} "
                sbatch_command += " -A {} --qos {} -p {} ".format(self.queue, self.queue, self.queue)

            # R1 case and R2 if needed
            READS = ["R1"]
            if self.is_paired(name) is True:
                READS += ["R2"]

            for RX in READS:
                print("  Merging {} ({} case)".format(name, RX))
                name.split("/")[-1]
                script_name = "{}_script.sh".format(os.path.split(name)[1])
                with open(script_name, "w") as fin:
                    cmd = self.get_pigz_cmd(name, RX)
                    fin.write(cmd)
                if which("sbatch") is not None:
                    cmd = sbatch_command.format(**{"thread": self.threads}) + " " + script_name
                else:
                    cmd = "sh {}".format(script_name)

                if dry_run is False:
                    process = subprocess.check_output(cmd.split())
                    processes.append(process)
                else:
                    processes.append("dryrun")

        if which("sbatch") is not None:
            if len(processes):
                print("Wait for those process to be over; type 'squeue | grep <your login>")
                for proc in processes:
                    print(
                        proc,
                    )

        self.processes = processes


@click.command(context_settings=CONTEXT_SETTINGS)
@click.option("--lanes", cls=OptionEatAll, required=True)
@click.option(
    "-o",
    "--output-directory",
    show_default=True,
    default="merging",
    required=True,
    help="Where to store the new fastq files",
)
@click.option(
    "--pattern",
    type=str,
    default="*/*fastq*gz",
    show_default=True,
    help="pattern for the input fastq files. Use quotes if wildcards are used",
)
@click.option("--threads", type=str, default=4, show_default=True, help="number of threads per job (pigz)")
@click.option("-s", "--use-sambamba", is_flag=True, help="use sambamba instad of samtools for sorting")
@click.option("--force", is_flag=True, default=False)
@click.option("--dry-run", is_flag=True, default=False)
@click.option("--slurm-queue", default=None)
def lane_merging(**kwargs):
    """Merge lanes

    This searches for data stored in this format:

              <sampleID_1>/*fastq.gz
              <sampleID_2>/*fastq.gz
              <sampleID_3>/*fastq.gz

          or::

              sampleID_L001_.fastq.gz
              sampleID_L002_.fastq.gz

    """
    lanes = [int(x) for x in eval(kwargs["lanes"])]
    print(lanes)

    c = LaneMerger(
        pattern=kwargs["pattern"],
        outdir=kwargs["output_directory"],
        queue=kwargs['slurm_queue'],
        lanes=lanes,
        force=kwargs["force"],
    )
    c.run(dry_run=kwargs["dry_run"])
