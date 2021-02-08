# -*- coding: utf-8 -*-
#
#  This file is part of Sequana software
#
#  Copyright (c) 2016 - Sequana Development Team
#
#  File author(s):
#    Thomas Cokelaer <thomas.cokelaer@pasteur.fr>
#
#  Distributed under the terms of the 3-clause BSD license.
#  The full license is in the LICENSE file, distributed with this software.
#
#  website: https://github.com/sequana/sequana
#  documentation: http://sequana.readthedocs.io
#
##############################################################################
"""Substract genomes from the raw reads"""
import os
import sys
import argparse
import glob
from subprocess import STDOUT
import subprocess

from easydev.console import purple

from sequana.scripts.tools import SequanaOptions
from sequana.bamtools import SAM
from sequana import FastQ
from sequana import logger

import colorlog
logger = colorlob.getLogger(__name__)


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
    def  __init__(self, prog="sequana_substractor"):
        usage = """%s reads (flag 256+4) saving the mapped reads in a file, and the unmapped in
another file\n""" % prog
        usage += """usage2: %s --input test.fastq --reference Phix174.fa\n""" % prog
        usage += """

        """
        super(Options, self).__init__(usage=usage, prog=prog,
                epilog=epilog,
                formatter_class=CustomFormatter)

        self.add_argument("--input", dest='input', type=str,
                            required=True, help="input FastQ file")
        self.add_argument("--output", dest='outfile', type=str,
                            default="reads.fastq", help="output FastQ filename")

        self.add_argument("--reference", dest="reference", type=str,
            default=None)
        self.add_argument("--references", dest="references", type=str,
            nargs="+", default=[])

        self.add_argument("--output-directory", dest='outdir', type=str,
                            default="sequana_substractor",
                            required=False, help="input fastq gzipped or not")
        self.add_argument("--mapper", dest='mapper', type=str,
                            default="minimap2", choices=["bwa", "minimap2"],
                            required=False, help="mapper minimap2 or bwa")

        self.add_threads(self)
        self.add_version(self)
        self.add_level(self)



class Substractor(object):
    def __init__(self, infile, references, outdir, mapper, threads=4):
        self.infile = infile
        self.references = references

        self.outdir = outdir
        self.threads = threads

        if os.path.exists(outdir):
            logger.info("using {} for output".format(outdir))
        else:
            os.mkdir(outdir)

        # this may be used later on for other mapper or methodology
        if mapper == "minimap2":
            self.mapper_cmd = "minimap2 -x map-pb -t {} {} {} -a > {}"
        elif mapper =="bwa" :
            self.mapper_cmd = "bwa mem -M -t {} {} {} > {}"
 
        f = FastQ(self.infile)
        self.L = len(f)
        logger.info("Found {} reads in input FastQ file\n\n".format(self.L))

    def run(self, output_filename):

        MAPPED = 0
        # temporary directory
        for i, reference in enumerate(self.references):
            if i == 0:
                infile = self.infile
            else:
                infile = outfile

            # we only accept reference ending in .fa or .fasta
            assert reference.endswith(".fa") or reference.endswith(".fasta")

            # keep only the basename
            outfile = os.path.basename(reference)
            outfile = outfile.replace(".fa", "").replace(".fasta", "")
            tag = outfile[0:8]
            outfile = "{}/mapping_{}.sam".format(self.outdir, tag)

            cmd = self.mapper_cmd.format(self.threads, reference, infile, outfile)

            # Now we need to extract the fastq from the SAM file.

            logger.info("Removing {}. Mapping starting".format(reference))
            logger.info(cmd)
            from subprocess import PIPE
            process = subprocess.call(cmd, shell=True, stderr=PIPE)

            results = self.splitter_mapped_unmapped(outfile, tag)

            # keep track of total mapped reads
            MAPPED += results["mapped"]

            outfile = "{}/{}.unmapped.fastq".format(self.outdir, tag)

            logger.info("{} mapped. {} reads remaining".format(
                results['mapped'], results["unmapped"]))
            print()

        # now we copy the last unmapped file into reads.fastq
        cmd = "cp {} {}".format(outfile, output_filename)
        process = subprocess.call(cmd, shell=True)

        logger.info("Your final file: {} with {} reads".format(
            output_filename, results['unmapped']))

        logger.info("all mapped and unmapped files: {}. Input was {}".format(
            MAPPED + results['unmapped'], self.L))

    def splitter_mapped_unmapped(self, filename, prefix):
        # helpful resources:
        # https://broadinstitute.github.io/picard/explain-flags.html
        logger.info("Creating 2 files (mapped and unmapped reads)")
        data = SAM(filename)

        results = {
            "flags": [],
            "mapped": 0,
            "unmapped": 0,
            "bad":0
        }
        logger.info("Please wait while creating output files")

        with open("{}/{}.unmapped.fastq".format(self.outdir, prefix), "w") as fnosirv:
            with open("{}/{}.mapped.fastq".format(self.outdir, prefix), "w") as fsirv:
                for a in data:
                    if a.flag & 2048: # suppl
                        # a bad read, we can just drop it
                        results['bad']+=1
                    elif a.flag & 1024: # PCR duplicate
                        results['bad']+=1
                    elif a.flag & 256: # secondary alignment
                        results["bad"] += 1
                    elif a.flag &16: # mapped
                        read = "@{}\n{}\n+\n{}\n".format(a.qname, a.query_sequence, a.qual)
                        assert len(a.query_sequence) == len(a.qual)
                        fsirv.write(read)
                        results["mapped"] += 1
                    elif a.flag & 4: # unmapped
                        read = "@{}\n{}\n+\n{}\n".format(a.qname, a.query_sequence, a.qual)
                        assert len(a.query_sequence) == len(a.qual)
                        fnosirv.write(read)
                        results["unmapped"] += 1
                    elif a.flag == 0: # mapped
                        read = "@{}\n{}\n+\n{}\n".format(a.qname, a.query_sequence, a.qual)
                        assert len(a.query_sequence) == len(a.qual)
                        fsirv.write(read)
                        results["mapped"] += 1
                    else:
                        logger.warning("{} flag not handled".format(a.flag))
                    results["flags"].append(a.flag)
        return results


def main(args=None):
    if args is None:
        args = sys.argv[:]

    print(purple("Welcome to sequana_substractor"))
    print(purple("WARNING. TESTED ON LONG READS ONLY. EXPERIMENTAL"))
    user_options = Options(prog="sequana_substractor")
    if len(args) ==1:
        args.append("--help")

    if "--version" in sys.argv:
        import sequana
        print(sequana.version)
        sys.exit(0)

    options = user_options.parse_args(args[1:])
    logger.setLevel(options.level)

    # build the references list
    references = []
    if options.reference:
        references.append(options.reference)
    if options.references:
        references = options.references
    options.references = references

    references = []
    # expand globs if any
    for ref in options.references:
        references.extend(glob.glob(ref))

    logger.info("{} references provided: {}".format(
        len(references), ",".join(references)))

    # call the entire machinery here
    sub = Substractor(options.input, references, options.outdir, 
        options.mapper, options.threads)
    sub.run(options.outfile)


if __name__ == "__main__":
    import sys
    main(sys.argv)

