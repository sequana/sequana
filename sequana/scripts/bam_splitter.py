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
"""Extract head of a zipped or unzipped FastQ file"""
import os
import sys
import argparse

from sequana.scripts.tools import SequanaOptions

from easydev.console import purple
from  collections import Counter


import colorlog
logger = colorlog.getLogger(__name__)




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
    def  __init__(self, prog="sequana_bam_splitter"):
        usage = """%s This BAM/SAM/CRAM filter removes unmapped and secondary
reads (flag 256+4) saving the mapped reads in a file, and the unmapped in
another file\n""" % prog
        usage += """usage2: %s --input yourbam.bam\n""" % prog
        usage += """usage2: %s --input yourbam.bam --prefix toto \n""" % prog
        usage += """usage2: %s --input yourbam.bam prefix --keep-unmapped\n""" % prog
        usage += """

        """
        super(Options, self).__init__(usage=usage, prog=prog,
                epilog=epilog,
                formatter_class=CustomFormatter)

        self.add_argument("--input", dest='input', type=str,
                            required=True, help="input SAM/BAM/CRAM file")
        self.add_argument("--output-directory", dest='outdir', type=str,
                            default=None,
                            required=False, help="input fastq gzipped or not")
        self.add_argument("--prefix", dest="prefix", default=None,
            required=False)
        self.add_argument("--keep-unmapped", dest="keep_unmapped",
                          action="store_true",
                          help="keep unmapped files")

        self.add_version(self)
        self.add_level(self)


def sniff(filename):
    logger.info("Sniffing file {}".format(filename))
    from sequana import BAM, SAM, CRAM
    from sequana.sniffer import sniffer
    datatype = sniffer(filename)
    if datatype == "SAM":
        logger.info("Input data in SAM format")
        data = SAM(filename)
    elif datatype == "BAM":
        logger.info("Input data in BAM format")
        data = BAM(filename)
    elif datatype == "CRAM":
        logger.info("Input data in CRAM format")
        data = CRAM(filename)
    else:
        raise ValueError("Your input file does not seem to be a valid SAM/BAM/CRAM file")
    return data


def splitter_mapped_unmapped(filename, prefix):
    logger.info("Creating 2 files (mapped and unmapped reads)")
    data = sniff(filename)

    count = 0
    flags = []
    match = 0
    unmatch = 0
    logger.info("Please wait while creating output files")
    with open("{}.unmapped.fastq".format(prefix), "w") as fnosirv:
        with open("{}.mapped.fastq".format(prefix), "w") as fsirv:
            for a in data:
                if a.flag & 256:
                    unmatch += 1
                elif a.flag & 4:
                    read = "@{}\n{}\n+\n{}\n".format(a.qname, a.query_sequence, a.qual)
                    assert len(a.query_sequence) == len(a.qual)
                    fnosirv.write(read)
                    unmatch += 1
                else:
                    read = "@{}\n{}\n+\n{}\n".format(a.qname, a.query_sequence, a.qual)
                    assert len(a.query_sequence) == len(a.qual)
                    fsirv.write(read)
                    match += 1
                flags.append(a.flag)
    return match, unmatch, flags


def splitter_mapped_only(filename, prefix):
    logger.info("Creating 1 file (mapped reads only). ")
    logger.info("Use --keep-unmapped to save unmapped reads.")
    data = sniff(filename)

    count = 0
    flags = []
    match = 0
    unmatch = 0

    logger.info("Please wait while creating output file")
    with open("{}.mapped.fastq".format(prefix), "w") as fsirv:
        for a in data:
            if a.flag & 260:
                unmatch += 1
                count += 1
            else:
                read = "@{}\n{}\n+\n{}\n".format(a.qname, a.query_sequence, a.qual)
                assert len(a.query_sequence) == len(a.qual)
                fsirv.write(read)
                match += 1
            flags.append(a.flag)
    return match, unmatch, flags


def _main(filename, prefix, keep_unmapped=True):
    if keep_unmapped:
        match, unmatch, flags = splitter_mapped_unmapped(filename, prefix)
    else:
        match, unmatch, flags = splitter_mapped_only(filename, prefix)
    return match, unmatch, flags


def main(args=None):
    if args is None:
        args = sys.argv[:]

    print(purple("Welcome to sequana_bam_splitter"))
    user_options = Options(prog="sequana_bam_splitter")
    if len(args) ==1:
        args.append("--help")


    if "--version" in sys.argv:
        import sequana
        print(sequana.version)
        sys.exit(0)

    options = user_options.parse_args(args[1:])

    # set the level
    logger.setLevel(options.level)
    logger.info("This SAM/BAM/CRAM splitter is used for paired or un-paired "
                "reads with perfectly mapped or unmapped reads (flags 0, 4, "
                "16). Others are dropped.")

    logger.info("Reading {}".format(options.input))

    # What prefix used for the output filename ?
    if options.prefix is None:
        prefix = options.input.rstrip(".bam")
        prefix = "test"
    else:
        prefix = options.prefix

    if options.outdir:
        prefix = options.outdir + os.sep + prefix
        if os.path.exists(options.outdir) is False:
            from easydev import mkdirs
            logger.info("Creating {} directory".format(options.outdir))
            mkdirs(options.outdir)


    match, unmatch, flags = _main(options.input, prefix,
        keep_unmapped=options.keep_unmapped)

    logger.info("Matched: {}".format(match))
    logger.info("Unmatched (flag 4 and 256): {}".format(unmatch))
    logger.info("All flags: {}".format(Counter(flags)))


if __name__ == "__main__":
    import sys
    main(sys.argv)

