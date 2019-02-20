# -*- coding: utf-8 -*-
#
#  This file is part of Sequana software
#
#  Copyright (c) 2016 - Sequana Development Team
#
#  File author(s):
#      Thomas Cokelaer <thomas.cokelaer@pasteur.fr>
#
#  Distributed under the terms of the 3-clause BSD license.
#  The full license is in the LICENSE file, distributed with this software.
#
#  website: https://github.com/sequana/sequana
#  documentation: http://sequana.readthedocs.io
#
##############################################################################
"""Merge Lanes. Be aware that this script works only on a SLURM cluster 

A future version will be more generic.

This only to be used by the biomics platform. You can of course contribute to
the improvment of this script

For biomics users, please see:
https://github.com/biomics-pasteur-fr/bioinfo/wiki/NextSeq:-how-to-merge-fusion-lanes
"""
import glob
import os
import subprocess
from subprocess import PIPE
import argparse



class Common():
    def __init__(self, pattern):
        """

        :param pattern: example */*fastq.gz

        """
        self.pattern = pattern
        self.discover_files()

    def discover_files(self, pattern=None):
        """Extract unique sample names"""
        if pattern:
            filenames = glob.glob(pattern)
        else:
            filenames = glob.glob(self.pattern)

        if len(filenames) == 0:
            ValueError("No files found. Are you in the correct path ? ")

        self.sampleIDs = sorted(list(set([x.split("_L00")[0] for x in filenames])))


class BackSpace(Common):
    """
    In backspace, the sample sheet contains the project, samplesID and name.

    The files are stored in a tree structure as follows::

        <Projet>/<sampleID>/<name>_R1_.fastq.gz
        <Projet>/<sampleID>/<name>_R2_.fastq.gz

    This script works at the sampleID level. You must be in the project
    directory.
    """
    def __init__(self, pattern="*/*fastq.gz", outdir="fusion", threads=4):
        super(BackSpace, self).__init__(pattern)
        self._outdir = None
        self.outdir = outdir
        self.threads = threads

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
            os.mkdir(self.outdir)

    def is_paired(self, name): 
        #filenames = glob.glob("{}_L*/{}*fastq.gz".format(name, name))    
        filenames = glob.glob("{}*/*fastq.gz".format(name))    
        R1 = sum([1 for filename in filenames if '_R1_' in filename])
        R2 = sum([1 for filename in filenames if '_R2_' in filename])

        if R1 == 4 and R2 == 4:
            return True
        elif R1 == 4 and R2 == 0:
            return False
        else:
            raise ValueError("Sample {} issue. Found {} R1 and {} R2. Expected 4 and 4 for paired data and 4 and 0 for single read".format(name, R1, R2))

    def get_pigz_cmd(self, sampleID, RX):
        assert RX in ['R1', 'R2']
        params = {"thread": self.threads, "sampleID": sampleID, "RX": RX}

        names = [x.split("/")[-1].split("_L00")[0] for x in 
                    glob.glob(sampleID + "*/*"+RX+"*")]
        print("  Found {} files for sample {} ({})".format(len(names), sampleID, RX))
        assert len(set(names)) == 1, "Non unique sample names {}".format(names)
        params['name'] = names[0]
        params['outdir'] = self.outdir

        output = "{outdir}/{sampleID}/{name}_{RX}_.fastq".format(**params)

        # Here, we should get 4 files with identical NAME (prefix), which should
        # be used for the output.
        #53_S53_L001_R2_001.fastq.gz

        cmd = "#!/bin/sh"
        cmd += "\nunpigz -p {thread} -c {sampleID}*/*_{RX}_*fastq.gz > " + output
        cmd += "\npigz --force -p {thread} " + output
        cmd  = cmd.format(**params)
        return cmd

    def run(self):
        import os
        print("Number of samples: {}".format(len(self.sampleIDs)))
        print(self.sampleIDs)
        print()

        processes = []
        for name in self.sampleIDs:
            print("{}, paired is {}".format(name, self.is_paired(name)))

            # create output directory
            import os
            if os.path.exists("{}/{}".format(self.outdir, name)):
                pass
            else:
                os.mkdir("{}/{}".format(self.outdir, name))

            # R1. Note the usage of wrap using --wrap " your command"
            sbatch_command = "sbatch -c {thread} --A biomics --qos biomics -p biomics"
            sbatch_command = "sbatch -c {thread} --qos fast"

            # R1 case and R2 if needed
            READS = ['R1']
            if self.is_paired(name) is True:
                READS += ['R2']
 
            for RX in READS:
                print("  fusionning {} ({} case)".format(name, RX))
                with open('script.sh', 'w') as fin:
                    cmd = self.get_pigz_cmd(name, RX)
                    fin.write(cmd)
                from shutil import which
                if which("sbatch") is not None:
                    cmd = sbatch_command.format(**{'thread': self.threads}) + ' script.sh'
                else:
                    cmd = "sh script.sh" 
                from subprocess import STDOUT
                process = subprocess.check_output(cmd.split())
                print(process)
                #proc = process.split()[-1]
                processes.append(process)

        if len(processes):
            print("Wait for those process to be over; type 'squeue | grep <your login>")
            for proc in processes:
                print(proc,)
        else:
            print('Found no fastq in ./* directories')
        self.processes = processes


class Options(argparse.ArgumentParser):
    def __init__(self, prog="backspace_fusion"):
        usage = """%s 

        backspace_fusion 

        This searches for data stored in this format:

            <sampleID_1>/*fastq.gz
            <sampleID_2>/*fastq.gz
            <sampleID_3>/*fastq.gz

        """.format(prog)
        super(Options, self).__init__(usage=usage, prog=prog,
                description="")

        # options to fill the config file
        self.add_argument("--output-directory", dest="outdir", type=str,
            default="fusion",
            help="Where to store the new fastq files")
        self.add_argument("--pattern", dest="pattern", type=str,
            default="*/*fastq*gz",
            help="Where to store the new fastq files")
        self.add_argument("--threads", dest="threads", type=str,
            defaut=4,
            help="number of threads per job (pigz)")
        #self.add_argument("--force", dest="force",
        #self.add_argument("--force", dest="force",
        #    default=False,
        #    action="store_true", help="""overwrite files""")


def main(args=None):

    user_options = Options(prog="sequana_backspace_fusion")
    if args is None:
        args = sys.argv

    # If --help or no options provided, show the help
    if "--help" in args:
        user_options.parse_args(["prog", "--help"])
    else:
        options = user_options.parse_args(args[1:])


    c = BackSpace(pattern=options.pattern, outdir=options.outdir)
    c.run()


if __name__ == "__main__":
    import sys
    main(sys.argv)
