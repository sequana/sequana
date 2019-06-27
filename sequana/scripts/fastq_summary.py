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
"""Extract head of a zipped or unzipped FastQ file"""
from sequana.fastq import FastQ
import glob
import sys
from optparse import OptionParser
import argparse


class Options(argparse.ArgumentParser):
    def  __init__(self, prog="sequana_fastq_summary"):
        usage = """%s --pattern "*fastq.gz" \n""" % prog
        usage += """Examples:

            sequana_fastq_summary --input "*R1*.fastq.gz" 

        """
        super(Options, self).__init__(usage=usage, prog=prog)
        self.add_argument("--pattern", dest='pattern', type=str,
                            required=True, help="input fastq gzipped or not")
 
def main(args=None):
    if args is None:
        args = sys.argv[:]

    user_options = Options(prog="sequana_fastq_count")

    if len(args) == 1 or "--help" in args:
        user_options.parse_args(["prog", "--help"])
    elif len(args) == 2:
        class SimpleOpt():
            pass
        options = SimpleOpt()
        options.input_filename = args[1]
    else:
        options = user_options.parse_args(args[1:])


    filenames = glob.glob(options.pattern)
    for filename in filenames:
        f = FastQ(filename)
        print(filename, len(f))


if __name__ == "__main__":
   import sys
   main(sys.argv)

