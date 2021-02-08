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
from sequana.vcf_filter import VCF

import sys
import argparse

from sequana.scripts.tools import SequanaOptions

import colorlog
logger = colorlog.getLogger(__name__)

from easydev.console import purple



class CustomFormatter(argparse.ArgumentDefaultsHelpFormatter,
                      argparse.RawDescriptionHelpFormatter):
    pass


epilog = purple("""
----

AUTHORS: Thomas Cokelaer
Documentation: http://sequana.readthedocs.io
Issues: http://github.com/sequana/sequana
        """)


class Options(argparse.ArgumentParser, SequanaOptions):
    def  __init__(self, prog="sequana_vcf_filter"):
        usage = """%s Only for VCF using mpileup version 4.1 for now\n""" % prog
        usage += """usage2: %s vcf_filter""" % prog
        usage += """Examples:

    sequana_vcf_filter --input test.vcf --quality 40
                --filter "AF1>0.95&AF1<0.05"
                --filter "MQ<30"


    You should look into the VCF file to figure out the valid TAGs. Then, you
    can apply various filters. 

    A filter should be interpreted as :

    ''filter out variants that agree with the filter''

    For example::

        --filter "DP<30"

    means ''remove SNPs with DP below 30''. We accept those types of comparison:

        DP<30
        DP<=30
        DP>30
        DP>=30

    For some tags, you want to keep values within or outside a range of
    values. You can then use the & and | characters::

        DP<30|>60  # to keep only values in the ranges [0-30] and [60-infinite]

    or

        DP>30&<60  # to keep only values in the range [30-60]

    Some tags stores a list of values. For instance DP4 contains 4 values.
    To filter the value at position 1, use e.g.::

        DP4[0]<0.5

    you can use the same convention for the range as above::

        DP4[0]>0.05&<0.95

    you may also need something like:

        sum(DP4[2]+DP4[3]) <2


    Note that you must use quotes to surround the filter values.

        """
        super(Options, self).__init__(usage=usage, prog=prog,
                epilog=epilog,
                formatter_class=CustomFormatter)

        self.add_argument("--input", dest='input_filename', type=str,
                            required=True, help="input fastq gzipped or not")

        self.add_argument("--quality", dest="quality",
                          type=int, default=-1, help="filter sites below this quality")

        self.add_argument("--apply-indel-filter", dest="apply_indel_filter",
                          action="store_true", help="remove INDELs")

        self.add_argument("--apply-dp4-filter", dest="apply_dp4_filter",
                          action="store_true", 
                          help="apply DP4 filters. see online doc, vcf_filter section")

        self.add_argument("--apply-af1-filter", dest="apply_af1_filter",
                          action="store_true", 
                          help="apply AF1 filters. see online doc, vcf_filter section")

        self.add_argument("--minimum-af1", dest="minimum_af1", 
            type=float, default=0.95, help="default to 0.95")
        self.add_argument("--minimum-ratio", dest="minimum_ratio", 
            type=float, default=0.75, help="default to 0.75")
        self.add_argument("--minimum-depth", dest="minimum_depth", 
            type=float, default=4, help="default to 4")
        self.add_argument("--minimum_depth-strand", dest="minimum_depth_strand", 
            type=float, default=2, help="default to 2")


        self.add_argument("--filter", dest="filter", action="append",
                        nargs=1, type=str, default=[],
                        help="Provide as many filters as you want. See example above ")

        self.add_argument("--output", dest="output_filename",
                            default="remaining.vcf", type=str)

        self.add_argument("--output-filtered", dest="output_filtered_filename",
                            default="filtered.vcf", type=str)

        self.add_argument('--level', dest="level",
            default="INFO", choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"])

        self.add_version(self)


def main(args=None):
    if args is None:
        args = sys.argv[:]


    print("Welcome to sequana_vcf_filter")
    user_options = Options(prog="sequana_vcf_filter")

    if "--version" in args:
        import sequana
        print(sequana.version)
        sys.exit()
    elif len(args) == 1 or "--help" in args:
        user_options.parse_args(["prog", "--help"])
    elif len(args) == 2:
        class SimpleOpt():
            pass
        options = SimpleOpt()
        options.input_filename = args[1]
    else:
        options = user_options.parse_args(args[1:])


    # set the level
    logger.setLevel(options.level)

    vcf = VCF(options.input_filename)
    try:
        vcf.vcf.filter_dict['QUAL'] =  options.quality
    except:
        vcf.vcf.filter_dict = {}
        vcf.vcf.filter_dict['QUAL'] =  options.quality

    vcf.vcf.apply_indel_filter = options.apply_indel_filter
    vcf.vcf.apply_dp4_filter = options.apply_dp4_filter
    vcf.vcf.apply_af1_filter = options.apply_af1_filter
    vcf.vcf.dp4_minimum_depth = options.minimum_depth
    vcf.vcf.dp4_minimum_depth_strand = options.minimum_depth_strand
    vcf.vcf.dp4_minimum_ratio = options.minimum_ratio
    vcf.vcf.minimum_af1 = options.minimum_af1
    vcf.vcf.filter_dict['INFO'] = {}
    vcf.vcf.filter_dict['QUAL'] =  options.quality


    for this in options.filter:
        this = this[0]
        signs = [">", "<", ">=", "<="]
        for sign in signs:
            if sign in this:
                key, value = this.split(sign, 1)
                key = key.strip()
                value = sign.strip() + value.strip()
                vcf.vcf.filter_dict['INFO'][key] = value
                break


    logger.info(vcf.vcf.filter_dict)

    res = vcf.vcf.filter_vcf(options.output_filename,
                       output_filtered=options.output_filtered_filename)

    print()
    #print(res)
    return res


if __name__ == "__main__":
   import sys
   main(sys.argv)

