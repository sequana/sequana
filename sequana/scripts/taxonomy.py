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
"""Standalone dedicated to taxonomic content (kraken)"""
import os
import sys
import argparse

from easydev import DevTools
from sequana.modules_report.kraken import KrakenModule
from sequana import KrakenPipeline, KrakenSequential
from sequana import KrakenDownload
from sequana import sequana_config_path as scfg
from sequana.taxonomy import Taxonomy
from sequana import sequana_config_path as cfg
from sequana.utils import config

import colorlog
logger = colorlog.getLogger(__name__)


class Options(argparse.ArgumentParser):
    def  __init__(self, prog="sequana_taxonomy"):
        usage = """Welcome to SEQUANA - Taxonomy standalone

        This standalone takes as input one or two (paired end) FastQ files.
        Using Kraken, Krona and some codecs from Sequana, it tries to identify
        the taxon/organism that match each reads found in the input files.

        Kraken requires an input database. We provide 3 databases but you can
        also use Sequana to help you building a customised one.

        We provide a DB called toydb. It contains only a few
        measles viruses. Its size is only 32Mb and should be used for testing
        and examples only.

        Another database is the so-called minikraken provided by Kraken's
        authors. It is about 4Gb and contains viruses and bacteria only.

        A third database is built with Sequana and is about 8Gb. It is
        stored on Synapse website and you will need an account (synapse.org).
        It contains about 22,000 species: viruses, bacteria, homo sapiens, fungi,

        Each DB can be downloaded using:

            sequana_taxonomy --download toydb

        Then, you need to use this kind of command:

            sequana_taxonomy --file1 R1.fastq --file2 R2.fastq
                --databases /home/user/.config/sequana/kraken_toydb
                --show-html --thread 4

AUTHORS: Thomas Cokelaer
Documentation: http://sequana.readthedocs.io
Issues: http://github.com/sequana/sequana
        """
        description = """DESCRIPTION:
        """
        super(Options, self).__init__(usage=usage, prog=prog,
                description=description)

        # options to fill the config file
        self.add_argument("--file1", dest="file1", type=str,
            help="""R1 fastq file (zipped) """)
        self.add_argument("--file2", dest="file2", type=str,
            help="""R2 fastq file (zipped) """)
        self.add_argument("--databases", dest="databases", type=str,
            nargs="+",
            help="""Path to a valid Kraken database(s). If you do not have any, use
                --download option. You may use several, in which case, an
                iterative taxonomy is performed as explained in online sequana
                documentation (see HierarchicalKRaken""")
        self.add_argument("--output-directory", dest="directory", type=str,
            help="""name of the output directory""", default="taxonomy")
        self.add_argument("--keep-temp-files", default=False,
            action="store_true", dest="keep_temp_files",
            help="keep temporary files (hierarchical case with several databases")
        self.add_argument("--thread", dest="thread", type=int,
            help="""number of threads to use (default 4)""", default=4)
        self.add_argument("--show-html", dest="html",
            action="store_true",
            help="""Results are stored in report/ directory and results are
                not shown by default""")
        self.add_argument("--download", dest="download", type=str,
            default=None, choices=["toydb"],
            help="""A toydb example to be downloaded.""")
        self.add_argument("--unclassified-out",
            help="save unclassified sequences to filename")
        self.add_argument("--classified-out",
            help="save classified sequences to filename")
        self.add_argument("--confidence", type=float, default=0, 
            help="confidence (kraken2 DB only)")
        self.add_argument("--update-taxonomy", action="store_true",
            default="False",
            help="Update the local NCBI taxonomy database to the last version")
        self.add_argument('--level', dest="level",
            default="INFO", choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"])


def main(args=None):

    if args is None:
        args = sys.argv[:]

    user_options = Options(prog="sequana")

    # If --help or no options provided, show the help
    if len(args) == 1:
        user_options.parse_args(["prog", "--help"])
    else:
       options = user_options.parse_args(args[1:])

    logger.setLevel(options.level)

    if options.update_taxonomy is True:
        tax = Taxonomy()
        logger.info("Will overwrite the local database taxonomy.dat in {}".format(cfg))
        tax.download_taxonomic_file(overwrite=True)
        sys.exit(0)

    # We put the import here to make the --help faster
    devtools = DevTools()

    if options.download:
        kd = KrakenDownload()
        kd.download(options.download)
        sys.exit()

    fastq = []
    if options.file1:
        devtools.check_exists(options.file1)
        fastq.append(options.file1)
    if options.file2:
        devtools.check_exists(options.file2)
        fastq.append(options.file2)

    if options.databases is None:
        logger.critical("You must provide a database")
        sys.exit(1)

    databases = []
    for database in options.databases:

        guess = os.sep.join([scfg, "kraken2_dbs", database])
        if os.path.exists(guess): # in Sequana path
            databases.append(guess)
        elif os.path.exists(database): # local database
            databases.append(database)
        else:
            msg = f"Invalid database name {database}. Neither found locally " + \
                  f"or in the sequana path {scfg}/kraken2_dbs; Use the --download option"
            logger.error(msg)
            raise IOError

    output_directory = options.directory + os.sep + "kraken"
    devtools.mkdirs(output_directory)

    # if there is only one database, use the pipeline else KrakenHierarchical
    _pathto = lambda x: f"{options.directory}/kraken/{x}" if x else x
    if len(databases) == 1:
        logger.info("Using 1 database")
        k = KrakenPipeline(fastq, databases[0], threads=options.thread,
            output_directory = output_directory, confidence=options.confidence)

        summary = k.run(output_filename_classified=_pathto(options.classified_out),
              output_filename_unclassified=_pathto(options.unclassified_out))
    else:
        logger.info("Using %s databases" % len(databases))
        k = KrakenSequential(fastq, databases,
            threads=options.thread,
            output_directory=output_directory + os.sep, 
            force=True,
            keep_temp_files=options.keep_temp_files,
            output_filename_unclassified=_pathto(options.unclassified_out),
            confidence=options.confidence)
        summary = k.run(output_prefix="kraken")

    with open(output_directory + "/summary.json", "w") as fh:
        import json
        json.dump( summary, fh, indent=4)

    # This statements sets the directory where HTML will be saved
    config.output_dir = options.directory

    # output_directory first argument: the directory where to find the data
    # output_filename is relative to the config.output_dir defined above
    kk = KrakenModule(output_directory, output_filename="summary.html")

    logger.info("Open ./%s/summary.html" % options.directory)
    logger.info("or ./%s/kraken/kraken.html" % options.directory)

    if options.html is True:
        ss.onweb()


if __name__ == "__main__":
   main(sys.argv)

