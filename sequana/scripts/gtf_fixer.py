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

import sys
from optparse import OptionParser
import argparse


class Options(argparse.ArgumentParser):
    def  __init__(self, prog="gtf_fixer"):
        usage = """%s input N output \n""" % prog
        usage += """usage2: %s gtf_filename""" % prog
        usage += """Examples:

            gtf_fixer --input test.gtf --output fixed.gtf

        """
        super(Options, self).__init__(usage=usage, prog=prog)
        self.add_argument("--input", dest='input', type=str,
                            required=True, help="input GTF file")
        self.add_argument("--output", dest='output', type=str,
                            required=True, help="output GTF file")

def main(args=None):
    if args is None:
        args = sys.argv[:]

    user_options = Options(prog="gtf_fixer")

    if len(args) == 1 or "--help" in args:
        user_options.parse_args(["prog", "--help"])
    elif len(args) == 2:
        class SimpleOpt():
            pass
        options = SimpleOpt()
        options.input_filename = args[1]
    else:
        options = user_options.parse_args(args[1:])

    from collections import defaultdict
    features = defaultdict(int)
    print("Scanning file")
    count = 0
    exon_ids = []
    with open(options.input, "r") as fin:
        with open(options.output, "w") as fout:
            line = fin.readline()
            count += 1
            while line:
                if line.startswith("#") or len(line.strip()) == 0:
                    fout.write(line)
                else:
                    entries = line.split(maxsplit=8)
                    features[entries[2]] += 1
                    try:
                        # sanity check that this is correctly written for all
                        # features
                        annotations = dict([x.strip().split(maxsplit=1) for x in entries[8].split(";") if len(x.strip())])
                    except:
                        print("warning line {}. could not parse annotations correctly".format(count))

                    # if exon feature, we want to make sure the ID are unique by
                    # appending the transcript_version if not already done
                    if entries[2] == "exon" and "." not in annotations["exon_id"]:
                        print("Fixing exon ID line {} for exon ID {}".format(count, annotations["exon_id"]))
                        exon_id = "{}.{}.{}".format(annotations['exon_id'].rstrip('"'), annotations['transcript_version'].strip('"'), annotations['transcript_id'].lstrip('"'))
                        newline = "\t".join(entries[0:8])+"\t"
                        # we do not resuse the dictionaru 'annotations' just to
                        # keep the order
                        for entry in entries[8].split(";"):
                            if len(entry.strip()):
                                x, y = entry.split(maxsplit=1)
                                if x.strip() != "exon_id":
                                    newline += "{} {};".format(x,y)
                                else:
                                    newline += "exon_id {};".format(exon_id)
                        newline +=";\n"
                        line = newline

                    fout.write(line)
                line = fin.readline()
                count += 1
                if count % 100000 == 0: print("scanned {} entries".format(count))

    print(features)

if __name__ == "__main__":
   import sys
   main(sys.argv)

