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
import os

from easydev import DevTools, execute, TempFile, md5

from sequana.lazy import pandas as pd
from sequana.lazy import pylab

from sequana.misc import wget
from sequana import sequana_config_path
from sequana import logger

logger.name = __name__


__all__ = ['KrakenResults', "KrakenPipeline", "KrakenAnalysis",
            "KrakenDownload", "KrakenHierarchical", "MultiKrakenResults"]


class KrakenResults(object):
    """Translate Kraken results into a Krona-compatible file


    If you run a kraken analysis with :class:`KrakenAnalysis`, you will end up
    with a file e.g. named kraken.out (by default).

    You could use kraken-translate but then you need extra parsing to convert
    into a Krona-compatible file. Here, we take the output from kraken and
    directly transform it to a krona-compatible file.

    ::

        k = KrakenResults("kraken.out")
        k.kraken_to_krona()

    Then format expected looks like::

        C    HISEQ:426:C5T65ACXX:5:2301:18719:16377    1    203    1:71 A:31 1:71
        C    HISEQ:426:C5T65ACXX:5:2301:21238:16397    1    202    1:71 A:31 1:71

    Where each row corresponds to one read.

    ::

        "562:13 561:4 A:31 0:1 562:3" would indicate that:

        the first 13 k-mers mapped to taxonomy ID #562
        the next 4 k-mers mapped to taxonomy ID #561
        the next 31 k-mers contained an ambiguous nucleotide
        the next k-mer was not in the database
        the last 3 k-mers mapped to taxonomy ID #562


    See kraken documentation for details.

    .. note:: a taxon of ID 1 (root) means that the read is classified but in
        differen domain. https://github.com/DerrickWood/kraken/issues/100

    .. note:: This takes care of fetching taxons and the corresponding lineages
        from online web services.

    """
    def __init__(self, filename="kraken.out"):
        """.. rubric:: **constructor**

        :param filename: the input from KrakenAnalysis class

        """
        self.filename = filename

        on_rtd = os.environ.get("READTHEDOCS", None) == "True"

        if on_rtd is False:
            from sequana.taxonomy import Taxonomy
            self.tax = Taxonomy(verbose=True)
            self.tax._load_flat_file() # make sure it is available locally
        else:
            class Taxonomy(object):
                from sequana import sequana_data # must be local
                df = pd.read_csv(sequana_data("test_taxon_rtd.csv"),
                        index_col=0)
                def get_lineage_and_rank(self, x):
                    # Note that we add the name as well here
                    ranks = ['kingdom', 'phylum', 'class', 'order',
                            'family', 'genus', 'species', 'name']
                    return [(self.df.loc[x][rank], rank) for rank in ranks]
            self.tax = Taxonomy()

        if filename:
            # This initialise the data
            self._parse_data()
            self._data_created = False

    def get_taxonomy_db(self, ids):
        """Retrieve taxons given a list of taxons

        :param list ids: list of taxons as strings or integers. Could also
            be a single string or a single integer
        :return: a dataframe

        .. note:: the first call first loads all taxons in memory and takes a
            few seconds but subsequent calls are much faster
        """
        # filter the lineage to keep only information from one of the main rank
        # that is superkingdom, kingdom, phylum, class, order, family, genus and
        # species
        ranks = ('kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species')

        if isinstance(ids, int):
            ids = [ids]

        if len(ids) == 0:
            return pd.DataFrame()

        logger.info('Retrieving taxon using sequana.Taxonomy')

        if isinstance(ids, list) is False:
            ids = [ids]

        lineage = [self.tax.get_lineage_and_rank(x) for x in ids]
        # Now, we filter each lineage to keep only relevant ranks
        # We drop the 'no rank' and create a dictionary
        # Not nice but works for now
        results = []
        for i, this in enumerate(lineage):
            default = dict.fromkeys(ranks, ' ')
            for entry in this:
                if entry[1] in ranks:
                    default[entry[1]] = entry[0]
                elif entry[1] == "superkingdom":
                    default["kingdom"] = entry[0]
            # Scientific name is the last entry tagged has no_rank  following
            # species TODO (check this assumption)
            # e.g. 351680 and 151529 have same 7 ranks so to differenatiate
            # them, the scientific name should be used.
            # By default, we will take the last one. If species or genus, we
            # repeat the term
            try:
                default['name'] = this[-1][0]
            except:
                default['name'] = "root (ambigous kingdom)"
            results.append(default)

        df = pd.DataFrame.from_records(results)
        df.index = ids
        df = df[list(ranks) + ['name']]
        df.index = df.index.astype(int)

        return df

    def _parse_data(self):
        taxonomy = {}

        logger.info("Reading kraken data from {}".format(self.filename))
        columns = ["status", "taxon", "length"]
        # we select only col 0,2,3 to save memoty, which is required on very
        # large files
        
        try:
            # each call to concat in the for loop below
            # will take time and increase with chunk position.
            # for 15M reads, this has a big cost. So chunksize set to 1M
            # is better than 1000 and still reasonable in memory
            reader = pd.read_csv(self.filename, sep="\t", header=None,
                               usecols=[0,2,3], chunksize=1000000)
        except pd.errors.EmptyDataError:
            logger.warning("Empty files. 100% unclassified ?")
            self.unclassified = "?" # size of the input data set
            self.classified = 0
            self._df = pd.DataFrame([], columns=columns)
            self._taxons = self._df.taxon
            return
        except pd.errors.ParserError:
            #raise NotImplementedError  # this section is for the case
            #    #only_classified_output when there is no found classified read
            raise NotImplementedError

        for chunk in reader:
            try:
                self._df
                self._df = pd.concat([self._df, chunk])
            except AttributeError:
                self._df = chunk

        self._df.columns = columns

        count = sum(self._df.taxon == 1)
        if count:
            logger.warning("Found %s taxons with root ID (1)" % count)

        # This gives the list of taxons as index and their amount
        # above, we select only columns 0, 2, 3  the column are still labelled
        # 0, 2, 3 in the df
        self._taxons = self._df.groupby("taxon").size()
        try:
            self._taxons.drop(0, inplace=True)
        except:
            pass # 0 may not be there
        self._taxons.sort_values(ascending=False, inplace=True)

        category = self.df.groupby("status").size()

        if 'C' in category.index:
            self.classified = category['C']
        else:
            self.classified = 0

        if 'U' in category.index:
            self.unclassified = category['U']
        else:
            self.unclassified = 0

    def _get_taxons(self):
        try:
            return self._taxons
        except:
            self._parse_data()
            return self._taxons
    taxons = property(_get_taxons)

    def _get_df(self):
        try:
            return self._df
        except:
            self._parse_data()
            return self._df
    df = property(_get_df)

    def _get_df_with_taxon(self, dbname):

        # line 14500
        # >gi|331784|gb|K01711.1|MEANPCG[331784] Measles virus (strain Edmonston), complete genome

        df = self.get_taxonomy_db([int(x) for x in self.taxons.index])
        df['count'] = self.taxons.values
        df.reset_index(inplace=True)
        newrow = len(df)
        df.loc[newrow] = "Unclassified"
        df.loc[newrow, 'count'] = self.unclassified
        df.loc[newrow, 'index'] = -1
        df.rename(columns={"index":"taxon"}, inplace=True)
        df["percentage"] = df["count"] / df["count"].sum() * 100

        # Now get back all annotations from the database itself.
        filename = dbname + os.sep + "annotations.csv"
        if os.path.exists(filename):
            annotations = pd.read_csv(filename)
            annotations.set_index("taxon", inplace=True)

            df2 = annotations.loc[df.taxon][['ena', 'gi', 'description']]
            # There are duplicates so let us keep the first one for now
            df2 = df2.reset_index().drop_duplicates(subset="taxon",
                keep="first").set_index("taxon")
            self.df2 = df2
            self.df1 = df.set_index("taxon")
            df = pd.merge(self.df1, df2, left_index=True, right_index=True)
            df.reset_index(inplace=True)
            starter = ['percentage', 'name', 'ena', 'taxon', "gi", 'count']
            df = df[starter + [x for x in df.columns if x not in starter and
                                x!="description"] +  ["description"]]

            df['gi'] = [int(x) for x in df['gi'].fillna(-1)]
            from easydev import precision
            df['percentage'] = [str(precision(x,2)) for x in df['percentage']]
        else:
            starter = ['taxon', 'count', 'percentage']
            df = df[starter + [x for x in df.columns if x not in starter]]

        df.sort_values(by="percentage", inplace=True, ascending=False)
        return df

    def kraken_to_csv(self, filename, dbname):
        df = self._get_df_with_taxon(dbname)
        df.to_csv(filename, index=False)
        return df

    def kraken_to_json(self, filename, dbname):
        df = self._get_df_with_taxon(dbname)
        df.to_json(filename)
        return df

    def kraken_to_krona(self, output_filename=None, mode=None, nofile=False):
        """

        :return: status: True is everything went fine otherwise False
        """
        if output_filename is None:
            output_filename = self.filename + ".summary"
        taxon_to_find = list(self.taxons.index)
        if len(taxon_to_find) == 0:
            logger.warning("No reads were identified. You will need a more complete database")
            self.output_filename = output_filename
            with open(output_filename, "w") as fout:
                fout.write("%s\t%s" % (self.unclassified, "Unclassified"))
            return False

        # classified reads as root  (1)
        """try:
            logger.warning("Removing taxon 1 (%s values) " % self.taxons.iloc[1])
            logger.info("Found %s taxons " % len(taxon_to_find))
            taxon_to_find.pop(taxon_to_find.index(1))
        except:
            pass
        """

        if len(taxon_to_find) == 0:
            return False

        if mode != "adapters":
            df = self.get_taxonomy_db(taxon_to_find)
            self.lineage = [";".join(this) for this in df[df.columns[0:-1]].values]
            self.scnames = list(df['name'].values)  # do we need a cast ?
        else:
            # Let us get the known adapters and their identifiers
            from sequana.adapters import AdapterDB
            adapters = AdapterDB()
            adapters.load_all()

            self.scnames = []

            for taxon in self.taxons.index:
                if str(taxon) in [1, "1"]:
                    self.scnames.append('unknown')
                    continue

                if str(taxon) not in list(adapters.df.identifier):
                    self.scnames.append('unknown')
                    continue

                self.scnames.append(adapters.get_name(taxon))
            self.lineage = ["Adapters;%s"% x for x in self.scnames]

            assert len(self.lineage) == len(self.taxons)
            assert len(self.scnames) == len(self.taxons)

        # Now save the file
        self.output_filename = output_filename
        with open(output_filename, "w") as fout:
            for i, this in enumerate(self.lineage):
                taxon = taxon_to_find[i]
                count = self.taxons.loc[taxon]
                line = str(count)+"\t"+"\t".join(this.split(';'))
                line += " " +self.scnames[i]
                fout.write(line+'\n')
            try:
                fout.write("%s\t%s" % (self.unclassified, "Unclassified"))
            except:
                pass #unclassified may not exists if all classified
        self._data_created = True
        return True

    def plot(self, kind="pie", cmap="tab20c", threshold=1, radius=0.9,
                textcolor="red", **kargs):
        """A simple non-interactive plot of taxons

        :return: None if no taxon were found and a dataframe otherwise

        A Krona Javascript output is also available in :meth:`kraken_to_krona`

        .. plot::
            :include-source:

            from sequana import KrakenResults, sequana_data
            test_file = sequana_data("test_kraken.out", "testing")
            k = KrakenResults(test_file)
            df = k.plot(kind='pie')

        .. seealso:: to generate the data see :class:`KrakenPipeline`
            or the standalone application **sequana_taxonomy**.


        .. todo:: For a future release, we could use this kind of plot 
            https://stackoverflow.com/questions/57720935/how-to-use-correct-cmap-colors-in-nested-pie-chart-in-matplotlib
        """
        if len(self._df) == 0:
            return

        if self._data_created == False:
            status = self.kraken_to_krona()

        if kind not in ['barh', 'pie']:
            logger.error('kind parameter: Only barh and pie are supported')
            return
        # This may have already been called but maybe not. This is not time
        # consuming, so we call it again here

        if len(self.taxons.index) == 0:
            return None

        df = self.get_taxonomy_db(list(self.taxons.index))

        # we add the unclassified only if needed
        if self.unclassified > 0:
            df.loc[-1] = ["Unclassified"] * 8

        data = self.taxons.copy()

        # we add the unclassified only if needed
        if self.unclassified > 0:
            data.loc[-1] = self.unclassified


        data = data/data.sum()*100
        assert threshold > 0 and threshold < 100

        # everything below the threshold (1) is gather together and summarised
        # into 'others'
        others = data[data<threshold].sum()

        data = data[data>=threshold]
        names = df.loc[data.index]['name']

        data.index = names.values

        if others > 0:
            data.loc['others'] = others

        try:
            data.sort_values(inplace=True)
        except:
            data.sort(inplace=True)

        pylab.figure(figsize=(10, 8))
        pylab.clf()
        if kind == "pie":
            ax = data.plot(kind=kind, cmap=cmap, autopct='%1.1f%%',
                radius=radius, **kargs)
            pylab.ylabel(" ")
            for text in ax.texts:
                #  large, x-small, small, None, x-large, medium, xx-small,
                #  smaller, xx-large, larger
                text.set_size("small")
                text.set_color(textcolor)
            for wedge in ax.patches:
                wedge.set_linewidth(1)
                wedge.set_edgecolor("k")
            self.ax = ax
        elif kind == "barh":
            ax = data.plot(kind=kind,  **kargs)
            pylab.xlabel(" percentage ")

        return data

    def to_js(self, output="krona.html", onweb=False):
        if self._data_created == False:
            status = self.kraken_to_krona()
        execute("ktImportText %s -o %s" % (self.output_filename, output))
        if onweb is True:
            import easydev
            easydev.onweb(output)

    def boxplot_classified_vs_read_length(self):
        """Show distribution of the read length grouped by classified or not"""
        self.df[["status", "length"]].groupby('status').boxplot()


class KrakenPipeline(object):
    """Used by the standalone application sequana_taxonomy

    This runs Kraken on a set of FastQ files, transform the results
    in a format compatible for Krona, and creates a Krona HTML report.

    ::

        from sequana import KrakenPipeline
        kt = KrakenPipeline(["R1.fastq.gz", "R2.fastq.gz"], database="krakendb")
        kt.run()
        kt.show()

    .. warning:: We do not provide Kraken database within sequana. You may
        either download a database from https://ccb.jhu.edu/software/kraken/
        or use this class to download a toy example that will
        be stored in e.g .config/sequana under Unix platforms.
        See :class:`KrakenDownload`.

    .. seealso:: We provide a standalone application of this class, which is
        called sequana_taxonomy and can be used within a command shell.

    """
    def __init__(self, fastq, database, threads=4, output_directory="kraken",
            dbname=None):
        """.. rubric:: Constructor

        :param fastq: either a fastq filename or a list of 2 fastq filenames
        :param database: the path to a valid Kraken database
        :param threads: number of threads to be used by Kraken
        :param output_directory: output filename of the Krona HTML page
        :param dbname:

        Description: internally, once Kraken has performed an analysis, reads
        are associated to a taxon (or not). We then find the correponding
        lineage and scientific names to be stored within a Krona formatted file.
        KtImportTex is then used to create the Krona page.

        """
        # Set and create output directory
        self._devtools = DevTools()
        self.output_directory = output_directory
        self._devtools.mkdir(output_directory)
        self.ka = KrakenAnalysis(fastq, database, threads)

        if dbname is None:
            self.dbname = os.path.basename(database)
        else:
            self.dbname = dbname

    def run(self, output_filename_classified=None,
                output_filename_unclassified=None,
                only_classified_output=False):
        """Run the analysis using Kraken and create the Krona output

        .. todo:: reuse the KrakenResults code to simplify this method.

        """
        # Run Kraken (KrakenAnalysis)
        kraken_results = self.output_directory + os.sep + "kraken.out"

        self.ka.run(
            output_filename=kraken_results,
            output_filename_unclassified=output_filename_unclassified,
            output_filename_classified=output_filename_classified,
            only_classified_output=only_classified_output
        )

        # Translate kraken output to a format understood by Krona and save png
        # image
        self.kr = KrakenResults(kraken_results)

        df = self.kr.plot(kind="pie")
        pylab.savefig(self.output_directory + os.sep + "kraken.png")

        prefix = self.output_directory + os.sep

        self.kr.kraken_to_json(prefix + "kraken.json", self.dbname)
        self.kr.kraken_to_csv(prefix + "kraken.csv", self.dbname)

        # Transform to Krona HTML
        from snakemake import shell
        kraken_html = self.output_directory + os.sep + "kraken.html"
        status = self.kr.kraken_to_krona(output_filename=prefix+"kraken.out.summary")
        if status is True:
            shell("ktImportText %s -o %s" % (prefix+"kraken.out.summary", kraken_html))
        else:
            shell("touch {}".format(kraken_html))

    def show(self):
        """Opens the filename defined in the constructor"""
        from easydev import onweb
        onweb(self.output)


class KrakenAnalysis(object):
    """Run kraken on a set of FastQ files

    In order to run a Kraken analysis, we firtst need a local database.
    We provide a Toy example. The ToyDB is downloadable as follows ( you will
    need to run the following code only once)::

        from sequana import KrakenDownload
        kd = KrakenDownload()
        kd.download_kraken_toydb()

    .. seealso:: :class:`KrakenDownload`  for more database and
        :class:`sequana.kraken_builder.KrakenBuilder` to build your own
        databases

    The path to the database is required to run the analysis. It has been
    stored in the directory ./config/sequana/kraken_toydb under Linux platforms
    The following code should be platform independent::

        import os
        from sequana import sequana_config_path
        database = sequana_config_path + os.sep + "kraken_toydb")

    Finally, we can run the analysis on the toy data set::

        from sequana import sequana_data
        data = sequana_data("Hm2_GTGAAA_L005_R1_001.fastq.gz", "data")
        ka = KrakenAnalysis(data, database=database)
        ka.run()

    This creates a file named *kraken.out*. It can be interpreted with
    :class:`KrakenResults`
    """
    def __init__(self, fastq, database, threads=4 ):
        """.. rubric:: Constructor

        :param fastq: either a fastq filename or a list of 2 fastq filenames
        :param database: the path to a valid Kraken database
        :param threads: number of threads to be used by Kraken
        :param output: output filename of the Krona HTML page

        :param return:

        """
        self._devtools = DevTools()
        self._devtools.check_exists(database)

        self.database = database
        self.threads = threads

        # Fastq input
        if isinstance(fastq, str):
            self.paired = False
            self.fastq = [fastq]
        elif isinstance(fastq, list):
            if len(fastq) == 2:
                self.paired = True
            else:
                self.paired = False
            self.fastq = fastq
        else:
            raise ValueError("Expected a fastq filename or list of 2 fastq filenames")

        for this in self.fastq:
            self._devtools.check_exists(database)

    def run(self, output_filename=None, output_filename_classified=None,
            output_filename_unclassified=None, only_classified_output=False):
        """Performs the kraken analysis

        :param str output_filename: if not provided, a temporary file is used
            and stored in :attr:`kraken_output`.
        :param str output_filename_classified: not compressed
        :param str output_filename_unclassified: not compressed

        """
        if output_filename is None:
            self.kraken_output = TempFile().name
        else:
            self.kraken_output = output_filename

        params = {
            "database": self.database,
            "thread": self.threads,
            "file1": self.fastq[0],
            "kraken_output": self.kraken_output,
            "output_filename_unclassified": output_filename_unclassified,
            "output_filename_classified": output_filename_classified,
            }

        if self.paired:
            params["file2"] = self.fastq[1]

        command = "kraken %(file1)s "

        if self.paired:
            command += " %(file2)s --paired"

        command += " -db %(database)s "
        command += " --threads %(thread)s --output %(kraken_output)s --out-fmt legacy"
        # for kraken <=1.0 --out-fmt did not exist
        #command += " --out-fmt legacy"

        if output_filename_unclassified:
            command +=  " --unclassified-out %(output_filename_unclassified)s "

        if only_classified_output is True:
            command += " --only-classified-output"

        if output_filename_classified:
            command +=  " --classified-out %(output_filename_classified)s "

        command = command % params
        # Somehow there is an error using easydev.execute with pigz
        from snakemake import shell
        shell(command)


class KrakenHierarchical(object):
    """Kraken Hiearchical Analysis

    This runs Kraken on a FastQ file with multiple k-mer databases in a
    hierarchical way. Unclassified sequences with the first database are input
    for the second, and so on.

    The input may be a single FastQ file or paired, gzipped or not. FastA are
    also accepted.


    """
    def __init__(self, filename_fastq, fof_databases, threads=1,
                 output_directory="./kraken_hierarchical/", 
                 keep_temp_files=False, output_filename_unclassified=None, 
                force=False):
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
            parameter to True.
        """
        # When running kraken in paired mode and saving the unclassified reads
        # in a file, the output file (fastq) contains both R1 and R2 so there
        # are concatenated in the same file. Actually, if there is R1 and R2,
        # there are concatenated as R1 N R2 (with the letter N as link). 
        # So, in the hiearchical search, paired case, the first iteration has 2
        # input files, must subsequent iterations will have only one file as
        # input, that is the output of the previous run (provided by
        # --unclassified-out option)
        self.filename_fastq = filename_fastq

        # input databases may be stored in a file
        if isinstance(fof_databases, str) and os.path.exists(fof_databases):
            with open(fof_databases, 'r') as fof:
                self.databases = [absolute_path.split('\n')[0] for absolute_path in fof.readlines()]
        # or simply provided as a list
        elif isinstance(fof_databases, list):
            self.databases = fof_databases[:]
        else:
            raise TypeError("input databases must be a list of valid kraken "
                            "databases or a file (see documebntation)")
        self.threads = threads
        self.output_directory = output_directory
        self.keep_temp_files = keep_temp_files

        # check if the output directory already exist
        try:
            os.mkdir(output_directory)
        except OSError:
            if os.path.isdir(output_directory) and force is False:
                logger.error('Output directory %s already exists' % output_directory)
                raise Exception
            elif force is True:
                logger.warning("Output directory %s already exists. You may "
                    "overwrite existing results" % output_directory)

        # list of input fastq files
        if isinstance(filename_fastq, list) and len(filename_fastq) in [1, 2]:
            self.inputs = filename_fastq[:]
        elif isinstance(filename_fastq, str):
            self.inputs = [filename_fastq]
        else:
            msg = "input file must be a string or list of 2 filenames"
            msg += "\nYou provided {}".format(filename_fastq)
            raise TypeError(msg)

        self.unclassified_output = output_filename_unclassified

    def _run_one_analysis(self, iteration):
        """ Run one analysis """
        db = self.databases[iteration]
        logger.info("Analysing data using database {}".format(db))

        # By default, the output contains only classified reads
        only_classified_output = True

        # a convenient alias
        _pathto = lambda x: self.output_directory + x

        # the output is saved in this file
        file_kraken_class = _pathto("kraken_%d.out" % iteration)
        output_filename_unclassified = _pathto("unclassified_%d.fastq" % iteration)
        file_fastq_unclass = _pathto("unclassified_%d.fastq" % iteration)



        if iteration == 0:
            inputs = self.inputs
        else:
            inputs = self._list_kraken_input[iteration-1]


        # if this is the last iteration (even if iteration is zero), save
        # classified and unclassified in the final kraken results.
        if iteration == len(self.databases) -1:
            only_classified_output = False


        analysis = KrakenAnalysis(inputs, db, self.threads)
        analysis.run(output_filename=file_kraken_class,
                     output_filename_unclassified=output_filename_unclassified,
                     only_classified_output=only_classified_output)

        self._list_kraken_input.append(file_fastq_unclass)
        self._list_kraken_output.append(file_kraken_class)

        if self.keep_temp_files:
            result = KrakenResults(file_kraken_class)
            result.to_js("%s/krona_%d.html" %(self.output_directory, iteration))

    def run(self, dbname="multiple", output_prefix="kraken_final"):
        """Run the hierachical analysis

        This method does not return anything but creates a set of files:

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
            status = self._run_one_analysis(iteration)
            last_unclassified = self._list_kraken_input[-1]
            stat = os.stat(last_unclassified)
            if stat.st_size == 0:
                break

        # concatenate all kraken output files
        file_output_final = self.output_directory + os.sep + "%s.out" % output_prefix
        with open(file_output_final, 'w') as outfile:
            for fname in self._list_kraken_output:
                with open(fname) as infile:
                    for line in infile:
                        outfile.write(line)

        # create html report
        logger.info("Analysing results")
        result = KrakenResults(file_output_final)

        # TODO: this looks similar to the code in KrakenPipeline. could be factorised
        result.to_js("%s%s%s.html" % (self.output_directory, os.sep, output_prefix))
        result.plot(kind="pie")
        pylab.savefig(self.output_directory + os.sep + "kraken.png")
        prefix = self.output_directory + os.sep
        result.kraken_to_json(prefix + "kraken.json", dbname)
        result.kraken_to_csv(prefix + "kraken.csv", dbname)

        # remove kraken intermediate files (including unclassified files)
        if self.unclassified_output:
            # Just cp the last unclassified file
            import shutil
            shutil.copy2(self._list_kraken_input[-1], self.unclassified_output)

        if not self.keep_temp_files:
            for f_temp in self._list_kraken_output:
                os.remove(f_temp)
            for f_temp in self._list_kraken_input:
                os.remove(f_temp)


class KrakenDownload(object):
    """Utility to download Kraken DB and place them in a local directory

    ::

        from sequana import KrakenDownload
        kd = KrakenDownload()
        kd.download('toydb')
        kd.download('minikraken')

    A large database (8Gb) is available on synapse and has the following DOI::

        doi:10.7303/syn6171000

    It can be downloaded manually or if you have a Synapse login
    (https://www.synapse.org), you can use::

        from sequana import KrakenDownload
        kd = KrakenDownload()
        kd.downloaded("sequana_db1")
    """
    dv = DevTools()
    def download(self, name, verbose=True):
        if name == "minikraken":
            self._download_minikraken(verbose=verbose)
        elif name == "toydb":
            self._download_kraken_toydb(verbose=verbose)
        elif name == "sequana_db1":
            self._download_sequana_db1(verbose=verbose)
        else:
            raise ValueError("name must be toydb or minikraken, or sequana_db1")


    def _download_kraken_toydb(self, verbose=True):
        """Download the kraken DB toy example from sequana_data into
        .config/sequana directory

        Checks the md5 checksums. About 32Mb of data
        """
        dv = DevTools()
        base = sequana_config_path + os.sep + "kraken_toydb"
        taxondir = base + os.sep + "taxonomy"
        dv.mkdir(base)
        dv.mkdir(taxondir)

        baseurl = "https://github.com/sequana/data/raw/master/"

        # download only if required
        logger.info("Downloading the database into %s" % base)

        md5sums = [
            "28661f8baf0514105b0c6957bec0fc6e",
            "97a39d44ed86cadea470352d6f69748d",
            "d91a0fcbbc0f4bbac918755b6400dea6",
            "c8bae69565af2170ece194925b5fdeb9"]
        filenames = [
            "database.idx",
            "database.kdb",
            "taxonomy/names.dmp",
            "taxonomy/nodes.dmp"]

        for filename, md5sum in zip(filenames, md5sums):
            url = baseurl + "kraken_toydb/%s" % filename
            filename = base + os.sep + filename
            if os.path.exists(filename) and md5(filename) == md5sum:
                logger.warning("%s already present" % filename)
            else:
                logger.info("Downloading %s" % url)
                wget(url, filename)

    def _download_minikraken(self, verbose=True):
        dv = DevTools()
        base = sequana_config_path + os.sep + ""
        taxondir = base + os.sep + "taxonomy"
        dv.mkdir(base)
        dv.mkdir(taxondir)

        logger.info("Downloading minikraken (4Gb)")

        filename = base + os.sep + "minikraken.tgz"
        if os.path.exists(filename) and md5(filename) == "30eab12118158d0b31718106785195e2":
            logger.warning("%s already present" % filename)
        else:
            wget("https://ccb.jhu.edu/software/kraken/dl/minikraken.tgz", filename)
        # unzipping. requires tar and gzip

    def _download_from_synapse(self, synid, target_dir):
        try:
            from synapseclient import Synapse
        except ImportError:
            raise ImportError("Please install synapseclient using 'pip install synapseclient'")
        try:
            self._synapse.get(synid, downloadLocation=target_dir)
        except:
            self._synapse = Synapse()
            self._synapse.login()
            self._synapse.get(synid, downloadLocation=target_dir)

    def _download_sequana_db1(self, verbose=True):
        dbname = "sequana_db1"
        from easydev import md5
        dir1 = sequana_config_path + os.sep + dbname
        dir2 = dir1 + os.sep + "taxonomy"
        self.dv.mkdir(dir1)
        self.dv.mkdir(dir2)

        logger.info("Downloading about 8Gb of data (if not already downloaded) from"
            " Synapse into %s" % dir1)

        from os.path import exists
        filename = dir1 + "ena_list.txt"
        if exists(filename) and md5(filename) == "a9cc6268f3338d1632c4712a412593f2":
            pass
        else:
            self._download_from_synapse('syn6171700', dir1)

        # database.idx
        filename = dir1 + "database.idx"
        if exists(filename) and md5(filename) == "2fa4a99a4f52f2f04c5a965adb1534ac":
            pass
        else:
            self._download_from_synapse('syn6171017', dir1)

        # database.kdb ; this one is large (8Gb)
        filename = dir1 + "database.kdb"
        if exists(filename) and md5(filename) == "ff698696bfc88fe83bc201937cd9cbdf":
            pass
        else:
            self._download_from_synapse('syn6171107', dir1)

        # Then, the taxonomy directory
        filename = dir1 + "names.dmp"
        if exists(filename) and md5(filename) == "10bc7a63c579de02112d125a51fd65d0":
            pass
        else:
            self._download_from_synapse('syn6171286', dir2)

        filename = dir1 + "nodes.dmp"
        if exists(filename) and md5(filename) == "a68af5a60434e2067c4a0a16df873980":
            pass
        else:
            self._download_from_synapse('syn6171289', dir2)

        filename = dir1 + "taxons.txt"
        if exists(filename) and md5(filename) == "e78fbb43b3b41cbf4511d6af16c0287f":
            pass
        else:
            self._download_from_synapse('syn6171290', dir2)
        logger.info('done. You should have a kraken DB in %s' % dir1)

        # The annotations
        wget("https://github.com/sequana/data/raw/master/sequana_db1/annotations.csv",
            dir1 + os.sep + "annotations.csv")



class MultiKrakenResults():

    def __init__(self, filenames, sample_names=None):
        self.filenames = filenames
        if sample_names is None:
            self.sample_names = list(range(1,len(filenames)+1))
        else:
            self.sample_names = sample_names

    def get_df(self):
        import pandas as pd
        data = {}
        for sample, filename in zip(self.sample_names, self.filenames):
            df = pd.read_csv(filename)
            df = df.groupby("kingdom")['percentage'].sum()
            # if a taxon is obsolete, the kingdom is empty. 
            # We will set the kingdom as Unclassified and raise a warning
            # if the count is > 5%
            if " " in df.index:
                percent = df.loc[" "]
                if percent > 5:
                    logger.warning("Found {}% of taxons in obsolete category".format(percent))
                if "Unclassified" in df.index:
                    df.loc['Unclassified'] += df.loc[' ']
                    df.drop(" ", inplace=True)
                else:
                    df.loc['Unclassified'] = df.loc[' ']
                    df.drop(" ", inplace=True)
            data[sample] = df

        df = pd.DataFrame(data)
        #df.to_json(output.data)
        df = df.sort_index(ascending=False)
        return df

    def plot_stacked_hist(self, output_filename=None, dpi=200, kind="barh", 
        fontsize=10, edgecolor="k", lw=1, width=1, ytick_fontsize=10):
        df = self.get_df()
        df.T.plot(kind=kind, stacked=True, edgecolor=edgecolor, lw=lw,
            width=width)
        ax = pylab.gca()
        positions = pylab.yticks()
        #ax.set_yticklabel(positions, labels, fontsize=ytick_fontsize)
        pylab.xlabel("Percentage (%)", fontsize=fontsize)
        pylab.ylabel("Sample index/name", fontsize=fontsize)
        pylab.yticks(fontsize=ytick_fontsize)
        pylab.legend(title="kingdom")
        pylab.xlim([0, 100])

        if output_filename:
            pylab.savefig(output_filename, dpi=dpi)





