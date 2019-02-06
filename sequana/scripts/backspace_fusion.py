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


def get_sample_names():
    """Extract unique sample names"""
    filenames = glob.glob("*/*fastq.gz")
    if len(filenames) == 0:
        ValueError("No files found. Are you in the correct path ? ")
    names = sorted(list(set([x.split("_L00")[0] for x in filenames])))
    return names


def is_paired(name): 
    filenames = glob.glob("{}_L*/{}*fastq.gz".format(name, name))    
    R1 = sum([1 for filename in filenames if '_R1_' in filename])
    R2 = sum([1 for filename in filenames if '_R2_' in filename])

    if R1 == 4 and R2 == 4:
        return True
    elif R1 == 4 and R2 == 0:
        return False
    else:
        raise ValueError("Sample {} issue. Found {} R1 and {} R2. Expected 4 and 4 for paired data and 4 and 0 for single read".format(name, R1, R2))


if os.path.exists("fusion"):
    pass
else:
    print("Creating ./fusion directory")
    os.mkdir("fusion")


def get_pigz_cmd(name, RX, thread=4):
    assert RX in ['R1', 'R2']
    params = {"thread": thread, "name": name, "RX":RX}
    output = "fusion/{name}/{name}_{RX}_.fastq".format(**params)
    cmd = "#!/bin/sh\n unpigz -p {thread} -c {name}*/*{name}*_{RX}_*fastq.gz > " + output
    cmd += "\npigz -p {thread} " + output
    cmd  = cmd.format(**params)
    return cmd


processes = []
for name in get_sample_names():
    print("{}, paired is {}".format(name, is_paired(name)))

    # create output directory
    import os
    if os.path.exists("fusion/{}".format(name)):
        pass
    else:
        os.mkdir("fusion/{}".format(name))


    # R1. Note the usage of wrap using --wrap " your command"
    sbatch_command = "sbatch -c {thread} --A biomics --qos biomics -p biomics"
    sbatch_command = "sbatch -c {thread} --qos fast"


    # R1 case and R2 if needed
    READS = ['R1']
    if is_paired(name) is True: 
        READS += ['R2']

    for RX in READS:
        print("fusionning {} ({} case)".format(name, RX))
        with open('script.sh', 'w') as fin:
            cmd = get_pigz_cmd(name, RX)
            fin.write(cmd)
        cmd = sbatch_command.format(**{'thread':4}) + ' script.sh'
        print(cmd)
        from subprocess import STDOUT
        process = subprocess.check_output(cmd.split())
        proc = process.split()[-1]
        processes.append(proc)

if len(processes):
    print("Wait for those process to be over; type 'squeue | grep <your login>")
    for proc in processes:
        print(proc,)

else:
    print('Found no fastq in ./* directories')
