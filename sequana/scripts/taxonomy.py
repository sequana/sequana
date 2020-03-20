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
import shutil
import glob
import sys
import argparse

from easydev import DevTools
from sequana import logger
from sequana.modules_report.kraken import KrakenModule

logger.name = __name__


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
            default=None, choices=["sequana_db1", "toydb", "minikraken"],
            help="""download an official sequana DB. The sequana_db1 is stored
                in a dedicated Synapse page (www.synapse.org). minikraken
                is downloaded from the kraken's author page, and toydb from
                sequana github.""")
        self.add_argument("--unclassified-out",
            help="save unclassified sequences to filename")
        self.add_argument("--classified-out",
            help="save unclassified sequences to filename")
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

    logger.level = options.level

    # We put the import here to make the --help faster
    from sequana import KrakenPipeline
    from sequana.kraken import KrakenHierarchical
    devtools = DevTools()

    if options.download:
        from sequana import KrakenDownload
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


    from sequana import sequana_config_path as scfg
    if options.databases is None:
        logger.critical("You must provide a database")
        sys.exit(1)

    databases = []
    for database in options.databases:
        if database == "toydb":
            database = "kraken_toydb"
        elif database == "minikraken":
            database = "minikraken_20141208"

        if os.path.exists(scfg + os.sep + database): # in Sequana path
            databases.append(scfg + os.sep + database)
        elif os.path.exists(database): # local database
            databases.append(database)
        else:
            msg = "Invalid database name (%s). Neither found locally "
            msg += "or in the sequana path %s; Use the --download option"
            raise ValueError(msg % (database, scfg))

    output_directory = options.directory + os.sep + "kraken"
    devtools.mkdirs(output_directory)

    # if there is only one database, use the pipeline else KrakenHierarchical
    _pathto = lambda x: "{}/kraken/{}".format(options.directory, x) if x else x
    if len(databases) == 1:
        logger.info("Using 1 database")
        k = KrakenPipeline(fastq, databases[0], threads=options.thread,
            output_directory=output_directory)

        k.run(output_filename_classified=_pathto(options.classified_out),
              output_filename_unclassified=_pathto(options.unclassified_out))
    else:
        logger.info("Using %s databases" % len(databases))
        k = KrakenHierarchical(fastq, databases,
            threads=options.thread,
            output_directory=output_directory + os.sep, 
            force=True,
            keep_temp_files=options.keep_temp_files,
            output_filename_unclassified=_pathto(options.unclassified_out))
        k.run(output_prefix="kraken")

    # This statements sets the directory where HTML will be saved
    from sequana.utils import config
    config.output_dir = options.directory

    # output_directory first argument: the directory where to find the data
    # output_filename is relative to the config.output_dir defined above
    kk = KrakenModule(output_directory, output_filename="summary.html")

    logger.info("Open ./%s/summary.html" % options.directory)
    logger.info("or ./%s/kraken/kraken.html" % options.directory)

    if options.html is True:
        ss.onweb()

if __name__ == "__main__":
   import sys
   main(sys.argv)

