#  This file is part of Sequana software
#
#  Copyright (c) 2016-2021 - Sequana Development Team
#
#
#  Distributed under the terms of the 3-clause BSD license.
#  The full license is in the LICENSE file, distributed with this software.
#
#  website: https://github.com/sequana/sequana
#  documentation: http://sequana.readthedocs.io
#
##############################################################################
import os
import shutil
import sys
from pathlib import Path, PosixPath

import colorlog
from colormap import Colormap
from easydev import TempFile, md5
from snakemake import shell

from sequana import sequana_config_path
from sequana.lazy import numpy as np
from sequana.lazy import pandas as pd
from sequana.lazy import pylab
from sequana.misc import wget
from sequana.kraken.analysis import KrakenDB, KrakenAnalysis, KrakenResults
logger = colorlog.getLogger(__name__)


__all__ = [
    "KrakenSequential",
]




class KrakenSequential(object):
    """Kraken Sequential Analysis

    This runs Kraken on a FastQ file with multiple k-mer databases in a
    sequencial way way. Unclassified sequences with the first database are input
    for the second, and so on.

    The input may be a single FastQ file or paired, gzipped or not. FastA are
    also accepted.


    """

    def __init__(
        self,
        filename_fastq,
        fof_databases,
        threads=1,
        output_directory="./kraken_sequential/",
        keep_temp_files=False,
        output_filename_unclassified=None,
        output_filename_classified=None,
        force=False,
        confidence=0,
    ):
        """.. rubric:: **constructor**

        :param filename_fastq: FastQ file to analyse
        :param fof_databases: file that contains a list of databases paths
            (one per line). The order is important. Note that you may also
            provide a list of datab ase paths.
        :param threads: number of threads to be used by Kraken
        :param output_directory: name of the output directory
        :param keep_temp_files: bool, if True, will keep intermediate files
            from each Kraken analysis, and save html report at each step
        :param bool force: if the output directory already exists, the
            instanciation fails so that the existing data is not overrwritten.
            If you wish to overwrite the existing directory, set this
            parameter to iTrue.
        """
        self.filename_fastq = filename_fastq
        self.confidence = confidence

        # input databases may be stored in a file
        if isinstance(fof_databases, str) and os.path.exists(fof_databases):
            with open(fof_databases, "r") as fof:
                self.databases = [absolute_path.split("\n")[0] for absolute_path in fof.readlines()]
        # or simply provided as a list
        elif isinstance(fof_databases, (list, tuple)):
            self.databases = fof_databases[:]
        else:
            raise TypeError(
                "input databases must be a list of valid kraken2 " "databases or a file (see documentation)"
            )

        self.databases = [KrakenDB(x) for x in self.databases]

        for d in self.databases:
            if d.version != "kraken2":
                logger.error(f"input database {d} is not valid kraken2 ")
                sys.exit(1)

        self.threads = threads
        self.output_directory = output_directory
        self.keep_temp_files = keep_temp_files

        # check if the output directory already exist
        try:
            os.mkdir(output_directory)
        except OSError:
            if os.path.isdir(output_directory) and force is False:
                logger.error("Output directory %s already exists" % output_directory)
                raise Exception
            elif force is True:
                logger.warning(
                    "Output directory %s already exists. You may " "overwrite existing results" % output_directory
                )

        # list of input fastq files
        if isinstance(filename_fastq, list) and len(filename_fastq) in [1, 2]:
            self.inputs = filename_fastq[:]
        elif isinstance(filename_fastq, str):
            self.inputs = [filename_fastq]
        else:
            msg = "input file must be a string or list of 2 filenames"
            msg += "\nYou provided {}".format(filename_fastq)
            raise TypeError(msg)

        if len(self.inputs) == 1:
            self.paired = False
        elif len(self.inputs) == 2:
            self.paired = True

        self.unclassified_output = output_filename_unclassified
        self.classified_output = output_filename_classified

    def _run_one_analysis(self, iteration):
        """Run one analysis"""
        db = self.databases[iteration]
        logger.info("Analysing data using database {}".format(db))

        # a convenient alias
        _pathto = lambda x: self.output_directory / x

        # the output is saved in this file
        if self.paired:
            # if paired, kraken2 expect a # and then will create 2 files (1 and 2
            # )
            # Note that kraken adds a _ before the # (1,2) so no need to add one
            output_filename_unclassified = _pathto("unclassified_%d#.fastq" % iteration)
            file_fastq_unclass = [
                _pathto("unclassified_%d_1.fastq" % iteration),
                _pathto("unclassified_%d_2.fastq" % iteration),
            ]
        else:
            output_filename_unclassified = _pathto("unclassified_%d.fastq" % iteration)
            file_fastq_unclass = _pathto("unclassified_%d.fastq" % iteration)

        if iteration == 0:
            inputs = self.inputs
        else:
            inputs = self._list_kraken_input[iteration - 1]

        # if this is the last iteration (even if iteration is zero), save
        # classified and unclassified in the final kraken results.
        if iteration == len(self.databases) - 1:
            only_classified_output = False
        else:
            only_classified_output = True

        file_kraken_out = self.output_directory / "kraken_{}.out".format(iteration)

        # The analysis itself
        analysis = KrakenAnalysis(inputs, db, self.threads, confidence=self.confidence)

        analysis.run(
            output_filename=file_kraken_out,
            output_filename_unclassified=output_filename_unclassified,
            only_classified_output=only_classified_output,
        )

        # save input/output files.
        self._list_kraken_input.append(file_fastq_unclass)
        self._list_kraken_output.append(file_kraken_out)

    def run(self, dbname="multiple", output_prefix="kraken_final"):
        """Run the sequential analysis

        :param dbname:
        :param output_prefix:
        :return: dictionary summarizing the databases names and
            classified/unclassied

        This method does not return anything creates a set of files:

        - kraken_final.out
        - krona_final.html
        - kraken.png  (pie plot of the classified/unclassified reads)

        .. note:: the databases are run in the order provided in the constructor.
        """
        # list of all output to merge at the end
        self._list_kraken_output = []
        self._list_kraken_input = []

        # Iteration over the databases
        for iteration in range(len(self.databases)):
            # The analysis itself
            status = self._run_one_analysis(iteration)

            last_unclassified = self._list_kraken_input[-1]

            # If everything was classified, we can stop here
            try:  # handle special case of emmpty FastQ file
                try:
                    stat = os.stat(last_unclassified)
                    if stat.st_size == 0:
                        break
                except TypeError:
                    stat = os.stat(last_unclassified[0])
                    if stat.st_size == 0:
                        break
            except FileNotFoundError:
                break

        # concatenate all kraken output files
        file_output_final = self.output_directory / f"{output_prefix}.out"
        with open(file_output_final, "w") as outfile:
            for fname in self._list_kraken_output:
                with open(fname) as infile:
                    for line in infile:
                        outfile.write(line)

        logger.info("Analysing final results")
        result = KrakenResults(file_output_final, verbose=False)

        try:
            result.histo_classified_vs_read_length()
            pylab.savefig(self.output_directory / "hist_read_length.png")
        except Exception as err:
            logger.warning("hist read length could not be computed")

        try:
            result.boxplot_classified_vs_read_length()
            pylab.savefig(self.output_directory / "boxplot_read_length.png")
        except Exception as err:
            logger.warning("hist read length could not be computed")

        # TODO: this looks similar to the code in KrakenPipeline. could be factorised
        result.to_js("%s.html" % (self.output_directory / output_prefix))
        try:
            result.plot2(kind="pie")
        except Exception as err:
            logger.warning(err)
            result.plot(kind="pie")
        pylab.savefig(self.output_directory / "kraken.png")
        prefix = self.output_directory
        result.kraken_to_json(prefix / "kraken.json", dbname)
        result.kraken_to_csv(prefix / "kraken.csv", dbname)

        # remove kraken intermediate files (including unclassified files)
        if self.unclassified_output:
            # Just cp the last unclassified file
            try:
                # single-end data (one file)
                shutil.copy2(self._list_kraken_input[-1], self.unclassified_output)
            except:
                for i, x in enumerate(self._list_kraken_input[-1]):
                    shutil.copy2(x, self.unclassified_output.replace("#", str(i + 1)))

        if self.classified_output:
            # Just cp the last classified file
            shutil.copy2(self._list_kraken_input[-1], self.classified_output)

        summary = {"databases": [x.name for x in self.databases]}
        total = 0
        classified = 0
        for f_temp, db in zip(self._list_kraken_output, self.databases):
            # In theory, the first N-1 DB returns only classified (C) read
            # and the last one contains both
            try:
                df = pd.read_csv(f_temp, sep="\t", header=None, usecols=[0])
                C = sum(df[0] == "C")
                U = sum(df[0] == "U")
            except pd.errors.EmptyDataError:
                # if no read classified,
                C = 0
                U = 0
            total += U
            total += C
            classified += C
            summary[db.name] = {"C": C}
            if U != 0:  # the last one
                summary["unclassified"] = U
        summary["total"] = total
        summary["classified"] = classified

        if not self.keep_temp_files:
            # kraken_0.out
            for f_temp in self._list_kraken_output:
                os.remove(f_temp)

            # unclassified
            for f_temp in self._list_kraken_input:
                if isinstance(f_temp, str):
                    try:
                        os.remove(f_temp)
                    except FileNotFoundError:
                        pass
                elif isinstance(f_temp, list):
                    for this in f_temp:
                        try:
                            os.remove(this)
                        except FileNotFoundError:
                            pass
        return summary


