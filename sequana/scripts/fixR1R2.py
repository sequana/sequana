# -*- coding: utf-8 -*-
#
#  This file is part of Sequana software
#
#  Copyright (c) 2019 - Sequana Development Team
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
"""Append /1 or /2 at the end of each reads in fastq files

Why: Some fastq created from BAM/SAM will loose the information about R1/R2
notion because the comments are lost. This tools append /1 at the end of the
name. Note that if there is a comment the /1 is added to the end of the comment,
not the name.

Sequana quality control used to create fastq without /1 /2 . This tools is a
hack to help people fixing this issue themselves with their own code. It is not
required anymore if you use sequana quality control pipeline since it was fixed
in the code. 

"""
import os, glob, sys


def compress(filename):
    print('compressing {}'.format(filename))
    cmd = "pigz {}".format(filename)
    retcode = os.system(cmd)
    assert retcode == 0


def decompress(filename):
    print('decompressing {}'.format(filename))
    cmd = "unpigz {}".format(filename)
    retcode = os.system(cmd)
    assert retcode == 0


def fix_identifier(zipped_filenames, directory):

    N = len(zipped_filenames)

    for i, filename in enumerate(zipped_filenames):

        print("--------------- {}/{} ".format(i+1, N))
        decompress(filename)

        filename2 = filename.replace(".gz", "")
        print("Fixing {}".format(filename2))

        # Fix it
        with open(filename2, "r") as fin:
            # use basename in case input is made of full pathnames.
            outfile = directory + "/" + os.path.basename(filename2)
            with open(outfile, "w") as fout:
                if "R1" in filename:
                    suffix = "/1"
                elif "R2" in filename:
                    suffix = "/2"
                else:
                    raise ValueError("Expecting R1 or R2 in the filename")

                for line in fin.readlines():
                    if line.startswith("@"):
                        line = line.strip() + "{}\n".format(suffix)
                    fout.write(line)

        compress(outfile)
        compress(filename2)

def main():

    if len(sys.argv) != 3:
        print( """usage: python fixR1R2.py input fix

- input: a pattern for input files ending in gz
- fix:  he directory where fixed files are copied""")
        sys.exit()

    directory = sys.argv[2]
    inputs = sys.argv[1]

    if os.path.exists(directory) is False: 
        os.mkdir(directory)

    zipped_filenames = glob.glob(inputs)

    fix_identifier(zipped_filenames, directory)


if __name__ == "__main__":
    main()







