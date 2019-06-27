# -*- coding: utf-8 -*-
#
#  This file is part of Sequana software
#
#  Copyright (c) 2016 - Sequana Development Team
#
#  File author(s):
#      Thomas Cokelaer <thomas.cokelaer@pasteur.fr>
#      Dimitri Desvillechabrol <dimitri.desvillechabrol@pasteur.fr>,
#          <d.desvillechabrol@gmail.com>
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

from optparse import OptionParser
import argparse
from sequana import logger

logger.name = __name__


class Options(argparse.ArgumentParser):
    def  __init__(self,  prog="fastq_head"):
        usage = """%s input --nlines 10000 output \n""" % prog
        usage += """usage2: %s --input input --nlines N --output output""" % prog
        usage += """Examples:

            fastq_head input.fastq.gz 10000 output.fastq.gz
            fastq_head input.fastq.gz 10000 output.fastq
            fastq_head input.fastq 10000 output.fastq.gz
            fastq_head input.fastq 10000 output.fastq

        you can also use named arguments::

            fastq_head --input input.fastq.gz --nlines 10000 --ouput output.fastq.gz

        """
        super(Options, self).__init__(usage=usage,  prog=prog)
        self.add_argument("--nlines", dest='N', type=int, required=True,
                          help="Number of lines to extract.")
        self.add_argument("--input", dest='input_filename', type=str,
                            required=True, help="input fastq gzipped or not")
        self.add_argument("--output", dest='output_filename', type=str,
                            required=True,
                            help="output file with .gz extension or not")

def main(args=None):
    import sys
    if args is None:
        args = sys.argv[:]

    user_options = Options(prog="fastq_head")
    if len(args) == 1:
        user_options.parse_args(["prog", "--help"])
    elif len(args) == 4:
        class SimpleOpt():
            pass
        options = SimpleOpt()
        options.input_filename = args[1]
        options.N = int(args[2])
        options.output_filename = args[3]
    else:
        options = user_options.parse_args(args[1:])

    f = FastQ(options.input_filename)
    f.extract_head(N=options.N,
                   output_filename=options.output_filename)


if __name__ == "__main__":
   import sys
   main(sys.argv)

