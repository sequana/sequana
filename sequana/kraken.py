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
import sys
import shutil

from easydev import DevTools, execute, TempFile, md5

from sequana.lazy import pandas as pd
from sequana.lazy import pylab

from sequana.misc import wget
from sequana import sequana_config_path
from sequana import logger

logger.name = __name__


__all__ = ['KrakenResults', "KrakenPipeline", "KrakenAnalysis",
            "KrakenDownload", "KrakenSequential", "MultiKrakenResults"]


class KrakenDB():
    """Class to handle a kraken DB"""
    def __init__(self, filename):

        if isinstance(filename, KrakenDB):
            filename = filename.path

        if os.path.exists(filename) is False:
            possible_path1 = sequana_config_path + "/" + filename
            possible_path2 = sequana_config_path + "/kraken2_dbs/" + filename
            possible_path3 = sequana_config_path + "/kraken1_dbs/" + filename
            if os.path.exists(possible_path1) is True:
                self.path = possible_path1
            elif os.path.exists(possible_path2) is True:
                self.path = possible_path2
            elif os.path.exists(possible_path3) is True:
                self.path = possible_path3
            else:
                msg = "{} not found locally or in {}. ".format(
                    filename, sequana_config_path)
                raise IOError(msg)
        else:
            self.path = os.path.abspath(filename)

        self.name = os.path.basename(self.path)

    def _get_database_version(self):
        if os.path.exists(self.path + os.sep + "hash.k2d"):
            return "kraken2"
        elif os.path.exists(self.path + os.sep + "database.kdb"):
            return "kraken1"
        else:
            logger.error("Looks like an invalid kraken database directory")
            return None
    version = property(_get_database_version)

    def __repr__(self):
        return self.name


class KrakenResults(object):
    """Translate Kraken results into a Krona-compatible file


    If you run a kraken analysis with :class:`KrakenAnalysis`, you will end up
    with a file e.g. named kraken.out (by default).

    You could use kraken-translate but then you need extra parsing to convert
    into a Krona-compatible file. Here, we take the output from kraken and
    directly transform it to a krona-compatible file.

    If you use sequantial analysis, can you use kraken-translate  or
    --use-names. what to do with the several outputs?

    kraken1 provides kraken-translate, kraken2 uses the --use-names that needs
    extra parsing as well.

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

    For kraken2, format is slighlty different since it depends on paired or not.
    If paired, ::

        C   read1 2697049 151|151 2697049:117 |:| 0:1 2697049:116

    See kraken documentation for details.

    .. note:: a taxon of ID 1 (root) means that the read is classified but in
        differen domain. https://github.com/DerrickWood/kraken/issues/100

    .. note:: This takes care of fetching taxons and the corresponding lineages
        from online web services.

    """
    def __init__(self, filename="kraken.out", verbose=True):
        """.. rubric:: **constructor**

        :param filename: the input from KrakenAnalysis class

        """
        self.filename = filename

        on_rtd = os.environ.get("READTHEDOCS", None) == "True"

        if on_rtd is False:
            from sequana.taxonomy import Taxonomy
            self.tax = Taxonomy(verbose=verbose)
            self.tax.download_taxonomic_file() # make sure it is available locally
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

        if isinstance(ids, list) is False:
            ids = [ids]

        lineage = [self.tax.get_lineage_and_rank(x) for x in ids]
        # Now, we filter each lineage to keep only relevant ranks
        # There are a few caveats though as explained hereafter

        # we merge the kingdom and superkingdom and subkingdom
        results = []
        for i, this in enumerate(lineage):
            default = dict.fromkeys(ranks, ' ')
            for entry in this:
                if entry[1] in ranks:
                    default[entry[1]] = entry[0]
            # if there is a superkingdom, overwrite the kingdom
            for entry in this:
                if entry[1] == "superkingdom":
                    default["kingdom"] = entry[0]
            if default["kingdom"] == " ":
                for entry in this:
                    if entry[1] == "subkingdom":
                        default["kingdom"] = entry[0]

            # in theory, we have now populated all ranks; 
            # Yet, there are several cases:
            # 1. all ranks are filled: perfect
            # 2. some ranks are empty: we fill them with a space.
            # 3. all ranks are empty: 
            #    a. this is the root
            #    b. this may be expected. e.g for an artifical sequence 
            #    c. all ranks below species are empty --> this is probably
            #       what we will get e.g. for plasmids

            # case 3.b
            if set([x[1] for x in this]) == {'no rank', 'species'}:
                # we can ignore the root and keep the others
                # if we end up with more than 6 entries, this is annoying
                # let us put a warning for now.
                count = 0
                for x in this:
                    if x[1] == "no rank" and x[0] != 'root':
                        default[ranks[count]] = x[0]
                        count += 1
                    if count >6:
                        logger.warning("too many no_rank in taxon{}".format(ids[i]))
                        break

            # for the name, we take the last entry, which is suppose to be the
            # scientific name found, so the scientific name of the taxon itself.
            # Note that this is not alwyas the species rank name
            default['name'] = this[-1][0]
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
        # we select only col 0,2,3 to save memory, which is required on very
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
        percentage = (count / len(self._df) * 100)
        if percentage >= 1:
            logger.warning(
                "Found {} taxons of classified reads with root ID (1) ({} %)".format(
                    count, round(percentage,2)))

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

        logger.debug(self.taxons.iloc[0:10])

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
        # This is not used anymore. let us put a warning for now
        filename = dbname + os.sep + "annotations.csv"
        if os.path.exists(filename):
            print("Using _get_df_with_taxon with valid annotation.csv in {}".format(filename))
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
        try:
            df.to_json(filename, indent=4, orient="records")
        except:
            df.to_json(filename, orient="records")
        return df

    def kraken_to_krona(self, output_filename=None, nofile=False):
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

        df = self.get_taxonomy_db(taxon_to_find)
        self.lineage = [";".join(this) for this in df[df.columns[0:-1]].values]
        self.scnames = list(df['name'].values)  # do we need a cast ?

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

    def plot2(self, kind="pie", fontsize=12):
        import numpy as np
        import matplotlib.pyplot as plt

        taxons = self.taxons.copy()
        if len(self.taxons.index) == 0:
            return None

        df = self.get_taxonomy_db(list(self.taxons.index))
        self.dd = df
        if self.unclassified > 0:
            df.loc[-1] = ['Unclassified'] * 8
            taxons[-1] = self.unclassified
        df['ratio'] = taxons / taxons.sum() * 100


        data_class = df.groupby(['kingdom', 'class']).sum()
        data_species = df.groupby(['kingdom', 'species']).sum()

        X = []
        Y = []
        Z = []
        labels = []
        zlabels, ztaxons = [], []
        kingdom_colors = []
        inner_colors = []
        inner_labels=[]
        species_colors = []
        taxons =  df['species'].reset_index().set_index('species')

        for kingdom in data_class.index.levels[0]:
            # kingdom info
            X.append(data_class.loc[kingdom].ratio.sum())

            # class info
            y = list(data_class.loc[kingdom].ratio.values)
            temp = data_class.loc[kingdom]
            y1 = temp.query('ratio>=0.5')
            y2 = temp.query('ratio<0.5')
            y = list(y1.ratio.values) + list(y2.ratio.values)
            inner_labels += list(y1.ratio.index) + [''] * len(y2.ratio)
            Y.extend(y)

            #species info
            temp = data_species.loc[kingdom]
            z1 = temp.query('ratio>=0.5')
            z2 = temp.query('ratio<0.5')
            z = list(z1.ratio.values) + list(z2.ratio.values)
            zlabels += list(z1.ratio.index) + [''] * len(z2.ratio)
            Z.extend(z)

            labels.append(kingdom)

            if kingdom == "Eukaryota":
                this_cmap = plt.cm.Purples
            elif kingdom == "Unclassified":
                this_cmap = plt.cm.Greys
            elif kingdom == "Bacteria":
                this_cmap = plt.cm.Reds
            elif kingdom == "Viruses":
                this_cmap = plt.cm.Greens
            elif kingdom == "Archea":
                this_cmap = plt.cm.Yellows
            else:
                this_cmap = plt.cm.Blues

            kingdom_colors.append(this_cmap(0.8))
            inner_colors.extend(this_cmap(np.linspace(.6,.2, len(y))))
            species_colors.extend(this_cmap(np.linspace(.6,.2, len(z))))


        fig, ax = pylab.subplots(figsize=(9.5,7))
        size = 0.2

        pct_distance=0
        w1, l1 = ax.pie(X, radius=1-2*size, colors=kingdom_colors,
           wedgeprops=dict(width=size, edgecolor='w'), labels=labels,
            labeldistance=0.4)

        w2, l2 = ax.pie(Y, radius=1-size, colors=inner_colors,
                labels=[x.replace('Unclassified', "") for x in inner_labels],
               wedgeprops=dict(width=size, edgecolor='w'),
            labeldistance=0.65)

        # labels can be long. Let us cut them 
        zlabels2 = []
        for this in zlabels:
            if len(this)>30:
                zlabels2.append(this[0:30] + "...") 
            else:
                zlabels2.append(this)

        w3, l3 = ax.pie(Z, radius=1, colors=species_colors,
               labels=[x.replace("Unclassified", "") for x in zlabels2],
               wedgeprops=dict(width=size, edgecolor='w'),
                labeldistance=0.9)

        ax.set(aspect="equal")
        pylab.subplots_adjust(right=1, left=0, bottom=0, top=1)
        pylab.legend(labels, title="kingdom", loc="upper right",
            fontsize=fontsize)
        import webbrowser
        mapper = {k:v for k,v in zip(zlabels, Z)}
        def on_pick(event):
            wedge = event.artist
            label = wedge.get_label()
            if mapper[label] >1:
                taxon = taxons.loc[label, "index"]
                webbrowser.open("https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id={}".format(taxon))
            else:
                wedge.set_color("white")
        for wedge in w3:
            wedge.set_picker(True)
        fig.canvas.mpl_connect("pick_event", on_pick)
        
        # this is used to check that everything was okay in the rules   
        return df


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
        self.dd = data
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

        # if paired and kraken2, there are | in length to separate both reads.
        # to simplify, if this is the case, we will just take the first read
        # length for now. 
        df = self.df.copy()
        try: # kraken2
            df.length = df.length.apply(lambda x: int(x.split("|")[0]))
        except:
            pass
        df[["status", "length"]].groupby('status').boxplot()
        return df


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
            dbname=None, confidence=0):
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
        self.database = database
        self.ka = KrakenAnalysis(fastq, database, threads, confidence=confidence)

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
        self.kr = KrakenResults(kraken_results, verbose=False)

        try:
            self.kr.plot2(kind="pie")
        except:
            self.kr.plot(kind="pie")
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

        database = KrakenDB(self.database)

        summary = {"database": [database.name]}
        summary[database.name] = {"C": int(self.kr.classified)}
        summary["U"] = int(self.kr.unclassified)
        summary["total"] = int(self.kr.unclassified + self.kr.classified)
        # redundant but useful and compatible with sequential approach
        summary["unclassified"] = int(self.kr.unclassified)
        summary["classified"] = int(self.kr.classified)

        return summary

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
    def __init__(self, fastq, database, threads=4, 
            confidence=0):
        """.. rubric:: Constructor

        :param fastq: either a fastq filename or a list of 2 fastq filenames
        :param database: the path to a valid Kraken database
        :param threads: number of threads to be used by Kraken
        :param confidence: for kraken2 only

        :param return:

        """
        self.database = KrakenDB(database)
        self.threads = threads
        self.confidence = confidence

        # Fastq input
        if isinstance(fastq, str):
            self.paired = False
            self.fastq = [fastq]
        elif isinstance(fastq, list):
            if len(fastq) == 2:
                self.paired = True
            elif len(fastq) == 1:
                self.paired = False
            else:
                raise IOError(("You must provide 1 or 2 files"))
            self.fastq = fastq
        else:
            raise ValueError("Expected a fastq filename or list of 2 fastq filenames")

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
            dirname = os.path.dirname(output_filename)
            if os.path.exists(dirname) is False:
                os.makedirs(dirname)

        # make sure the required output directories exit:
        # and that the output filenames ends in .fastq
        if output_filename_classified:
            assert output_filename_classified.endswith(".fastq")
            dirname = os.path.dirname(output_filename_classified)
            if os.path.exists(dirname) is False:
                os.makedirs(dirname)

        if output_filename_unclassified:
            assert output_filename_unclassified.endswith(".fastq")
            dirname = os.path.dirname(output_filename_unclassified)
            if os.path.exists(dirname) is False:
                os.makedirs(dirname)
 
        params = {
            "database": self.database.path,
            "thread": self.threads,
            "file1": self.fastq[0],
            "kraken_output": self.kraken_output,
            "output_filename_unclassified": output_filename_unclassified,
            "output_filename_classified": output_filename_classified,
            }

        if self.paired:
            params["file2"] = self.fastq[1]

        if self.database.version == "kraken2":
            logger.info("Using kraken2 ")
            command = "kraken2 --confidence {}".format(self.confidence)
        elif self.database.version == "kraken1":
            # for kraken <=1.0 --out-fmt did not exist
            #command += " --out-fmt legacy"
            # version 0.8.3 of sequana, we changed the output from legacy to
            # interleaved. legacy puts R1 and R2 together (separated by a N)
            # but N's may be present in the reads itself, so we prefer to use
            # the interleaved now so that it is easy to split the interleaved
            # R1/R2 into two files R1 and R2
            # we also added the --fastq-output to be compatible/coherent with kraken2
            # behaviour
            logger.info("Using kraken1. You should use version 1.1 ")
            command = "kraken "
            # --fastq-output seems to work only if paired
            # therefore --fastq-output can only be used with --paired and
            # --out-fmt interleaved makes sense only with --paired and therefore
            # --fastq-output. Consequently, single end data output will only be
            # available as fasta
            if self.paired:
                command += " --fastq-output --out-fmt interleaved "

        command += " %(file1)s"

        if self.paired:
            command += " %(file2)s --paired"

        command += " --db %(database)s "
        command += " --threads %(thread)s --output %(kraken_output)s "

        # If N is number of reads unclassified 3 cases depending on out-fmt
        # choice
        # case1 --paired and out-fmt legacy saved fasta R1 and R2 together on N lines
        # case2 --paired and out-fmt interleaved saved fasta R1 and R2 alternatively on 2N lines
        # case3 --paired and out-fmt paired saved R1 on N lines. Where is R2 ????
        # Note, that there is always one single file. So, the only way for
        # kraken to know that this new files (used as input) is paired, is to
        # use --paired.

        # In any case, this new file looks like an R1-only file. Indeed, if
        # interleaved, all data inside the file, if legacy, The R1 and R2 are
        # separated by N but a unique sequence. If --out-fmt is paired, this is
        # annoying. Indeed, half of the data is lost.

        # So, if now input is
        # case1, we cannot provide --paired
        # case2 we cannot either, so how are R1 and R2 taken care of ?
        # besides, if provided, the interleaved input is seen as single ended.
        # Indeed, if provided, --out-fmt cannot be interleaved since krakne1
        # complains that input is not paired.
        # case3, only R1 so we cannot use --paired

        # if kraken2, there is no --out-fmt option, so output is always a fastq
        # with either R1 only or two output files.
        # If we omit the --paired options, the 2 input R1 and R2 are considered
        # as 2 different unrelated samples 
        # if we use --paired we now must have # in the file name, and then
        # the two files are created
        if self.database.version == "kraken2":
            if output_filename_unclassified:
                command +=  " --unclassified-out %(output_filename_unclassified)s "
            if output_filename_classified:
                command +=  " --classified-out %(output_filename_classified)s "

        # For kraken1, if paired, we first save into a temporary file as interleaved data
        # and then split the data. That way, it is compatible with
        # kraken2 if databases are mixed. If not paired, nothing special to do;
        # just save the data
        if self.database.version == "kraken1":
            if output_filename_unclassified:
                command +=  " --unclassified-out %(output_filename_unclassified)s "
            if output_filename_classified:
                command +=  " --classified-out %(output_filename_classified)s "

        # This is for the kraken output, not classified fastq
        if only_classified_output is True:
            # This was a convenient option in kraken1, not available anymore in
            # kraken2. See later after the command call 
            command += " --only-classified-output"

        command = command % params
        # Somehow there is an error using easydev.execute with pigz
        from snakemake import shell
        logger.debug(command)
        shell(command)

        if self.database.version == "kraken2" and only_classified_output:
            # kraken2 has no classified_output option. we mimic it here below
            # just to get a temporary filename
            fout = TempFile()
            outname = fout.name
            newfile = open(outname, "w")
            with open(output_filename, "r") as fin:
                for line in fin.readlines():
                    if line.startswith("C"):
                        newfile.write(line)
            newfile.close()
            shutil.move(outname, output_filename)
            # now let us delete the temp file
            # The temp file has been moved so no need to delete it
            # fout.delete()

        # a simple utility function
        try:
            from itertools import izip_longest
        except:
            from itertools import zip_longest as izip_longest 
        def grouper(iterable):
            args = [iter(iterable)] * 8
            return izip_longest(*args)

        if self.database.version == "kraken1" and self.paired:
            logger.info("This is kraken1 paired, let us split the interleaved data set")
            if output_filename_unclassified:
                f1 = output_filename_unclassified.replace(".fastq", "_1.fastq")
                f2 = output_filename_unclassified.replace(".fastq", "_2.fastq")
                with open(f1, "w") as R1, open(f2, "w") as R2:
                    with open(output_filename_unclassified, "r") as fin:
                        for lines in grouper(fin):
                            R1.write("".join(lines[0:4]))
                            R2.write("".join(lines[4:]))
                # now we can remove the temporary file
                os.remove(output_filename_unclassified)
            if output_filename_classified:
                f1 = output_filename_classified.replace(".fastq", "_1.fastq")
                f2 = output_filename_classified.replace(".fastq", "_2.fastq")
                with open(f1, "w") as R1, open(f2, "w") as R2:
                    with open(output_filename_classified, "r") as fin:
                        for lines in grouper(fin):
                            R1.write("".join(lines[0:4]))
                            R2.write("".join(lines[4:]))
                # now we can remove the temporary file
                os.remove(output_filename_classified)


class KrakenSequential(object):
    """Kraken Sequential Analysis

    This runs Kraken on a FastQ file with multiple k-mer databases in a
    sequencial way way. Unclassified sequences with the first database are input
    for the second, and so on.

    The input may be a single FastQ file or paired, gzipped or not. FastA are
    also accepted.


    """
    def __init__(self, filename_fastq, fof_databases, threads=1,
                 output_directory="./kraken_sequential/",
                 keep_temp_files=False, 
                 output_filename_unclassified=None,
                 output_filename_classified=None,
                 force=False, confidence=0):
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
            with open(fof_databases, 'r') as fof:
                self.databases = [absolute_path.split('\n')[0] for absolute_path in fof.readlines()]
        # or simply provided as a list
        elif isinstance(fof_databases, list):
            self.databases = fof_databases[:]
        else:
            raise TypeError("input databases must be a list of valid kraken "
                            "databases or a file (see documebntation)")

        self.databases = [KrakenDB(x) for x in self.databases]

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

        if len(self.inputs) == 1:
            self.paired = False
        elif len(self.inputs) == 2:
            self.paired = True

        self.unclassified_output = output_filename_unclassified
        self.classified_output = output_filename_classified

    def _run_one_analysis(self, iteration):
        """ Run one analysis """
        db = self.databases[iteration]
        logger.info("Analysing data using database {}".format(db))


        # a convenient alias
        _pathto = lambda x: self.output_directory + x

        # the output is saved in this file
        if db.version == "kraken2" and self.paired:
            # if paired, kraken2 expect a # and then will create 2 files (1 and 2
            # )
            # Note that kraken adds a _ before the # (1,2) so no need to add one
            output_filename_unclassified = _pathto("unclassified_%d#.fastq" % iteration)
            file_fastq_unclass = [
                _pathto("unclassified_%d_1.fastq" % iteration),
                _pathto("unclassified_%d_2.fastq" % iteration)]
        elif db.version == "kraken2" and self.paired is False:
            output_filename_unclassified = _pathto("unclassified_%d.fastq" % iteration)
            file_fastq_unclass = _pathto("unclassified_%d.fastq" % iteration)
        elif db.version == "kraken1" and self.paired is False:
            output_filename_unclassified = _pathto("unclassified_%d.fastq" % iteration)
            file_fastq_unclass =  _pathto("unclassified_%d.fastq" % iteration)
        elif db.version == "kraken1" and self.paired is True:
            # if paired, kraken1 expects 1 file and creates 1 file.
            # To be compatible with inputs of kraken2, if paired, the output
            # must be made of two files that will be created by KrakenAnalysis
            # class by us (using the single file created by kraken1)
            output_filename_unclassified = _pathto("unclassified_%d.fastq" % iteration)
            file_fastq_unclass = [
                _pathto("unclassified_%d_1.fastq" % iteration),
                _pathto("unclassified_%d_2.fastq" % iteration)]

        if iteration == 0:
            inputs = self.inputs
        else:
            inputs = self._list_kraken_input[iteration-1]

        # if this is the last iteration (even if iteration is zero), save
        # classified and unclassified in the final kraken results.
        if iteration == len(self.databases) -1:
            only_classified_output = False
        else:
            only_classified_output = True

        file_kraken_out = self.output_directory + "/kraken_{}.out".format(iteration)

        # The analysis itself
        analysis = KrakenAnalysis(inputs, db, self.threads,
            confidence=self.confidence)

        analysis.run(output_filename=file_kraken_out,
                     output_filename_unclassified=output_filename_unclassified,
                     only_classified_output=only_classified_output)

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
        import shutil
        # list of all output to merge at the end
        self._list_kraken_output = []
        self._list_kraken_input = []

        # Iteration over the databases
        for iteration in range(len(self.databases)):
            #kraken_out = self.output_directory + "/kraken_{}.out".format(iteration)

            # The analysis itself
            status = self._run_one_analysis(iteration)

            # We save intermediate krona except last one that will be done
            # anyway
            if self.keep_temp_files and iteration != len(self.databases) - 1:
                kraken_out = self._list_kraken_output[-1] 
                result = KrakenResults(kraken_out, verbose=False)
                result.to_js("%s/krona_%d.html" %(self.output_directory, iteration))

            last_unclassified = self._list_kraken_input[-1]

            # If everything was classified, we can stop here
            if isinstance(last_unclassified, str):
                stat = os.stat(last_unclassified)
                if stat.st_size == 0:
                    break
            elif isinstance(last_unclassified, list):
                stat = os.stat(last_unclassified[0])
                if stat.st_size == 0:
                    break

        # concatenate all kraken output files
        file_output_final = self.output_directory + os.sep + "%s.out" % output_prefix
        with open(file_output_final, 'w') as outfile:
            for fname in self._list_kraken_output:
                with open(fname) as infile:
                    for line in infile:
                        outfile.write(line)

        logger.info("Analysing final results")
        result = KrakenResults(file_output_final, verbose=False)

        # TODO: this looks similar to the code in KrakenPipeline. could be factorised
        result.to_js("%s%s%s.html" % (self.output_directory, os.sep, output_prefix))
        try:
            result.plot2(kind="pie")
        except:
            result.plot(kind="pie")
        pylab.savefig(self.output_directory + os.sep + "kraken.png")
        prefix = self.output_directory + os.sep
        result.kraken_to_json(prefix + "kraken.json", dbname)
        result.kraken_to_csv(prefix + "kraken.csv", dbname)

        # remove kraken intermediate files (including unclassified files)
        if self.unclassified_output:
            # Just cp the last unclassified file
            shutil.copy2(self._list_kraken_input[-1], self.unclassified_output)

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
                df = pd.read_csv(f_temp, sep="\t", header=None,usecols=[0])
                C = sum(df[0] == "C")
                U = sum(df[0] == "U")
            except pd.errors.EmptyDataError:
                # if no read classified, 
                C = 0
                U = 0
            total+=U
            total += C
            classified += C
            summary[db.name] = {"C":C}
            if U != 0: # the last one
                summary['unclassified'] = U
        summary['total'] = total
        summary['classified'] = classified

        if not self.keep_temp_files:
            # kraken_0.out
            for f_temp in self._list_kraken_output:
                os.remove(f_temp)

            # unclassified
            for f_temp in self._list_kraken_input:
                if isinstance(f_temp, str):
                    os.remove(f_temp)
                elif isinstance(f_temp, list):
                    for this in f_temp:
                        os.remove(this)
        return summary
        

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
        elif name == "toydb2":
            self._download_kraken2_toydb(verbose=verbose)
        elif name == "sequana_db1":
            self._download_sequana_db1(verbose=verbose)
        else:
            raise ValueError("name must be toydb or minikraken, or sequana_db1")

    def _download_kraken2_toydb(self, verbose=True):
        """Download the kraken DB toy example from sequana_data into
        .config/sequana directory

        Checks the md5 checksums. About 32Mb of data
        """
        dv = DevTools()
        base = sequana_config_path + os.sep + "kraken2_dbs/toydb_nomasking"
        dv.mkdir(base)

        baseurl = "https://github.com/sequana/data/raw/master/"

        # download only if required
        logger.info("Downloading the database into %s" % base)

        md5sums = [
            "31f4b20f9e5c6beb9e1444805264a6e5",
            "733f7587f9c0c7339666d5906ec6fcd3",
            "7bb56a0f035b27839fb5c18590b79263"]

        filenames = [
            "hash.k2d",
            "opts.k2d",
            "taxo.k2d"]

        for filename, md5sum in zip(filenames, md5sums):
            url = baseurl + "kraken2_toydb/%s" % filename
            filename = base + os.sep + filename
            if os.path.exists(filename) and md5(filename) == md5sum:
                logger.warning("%s already present" % filename)
            else:
                logger.info("Downloading %s" % url)
                wget(url, filename)

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
    """Select several kraken output and creates summary plots

    ::

        import glob
        mkr = MultiKrakenResults(glob.glob("*/*/kraken.csv"))
        mkr.plot_stacked_hist()

    """
    def __init__(self, filenames, sample_names=None):
        self.filenames = filenames
        if sample_names is None:
            self.sample_names = list(range(1,len(filenames)+1))
        else:
            self.sample_names = sample_names

    def get_df(self, limit=5):
        import pandas as pd
        data = {}
        for sample, filename in zip(self.sample_names, self.filenames):
            df = pd.read_csv(filename)
            if "kingdom" not in df.columns:
                for name in "kingdom,phylum,class,order,family,genus,species,name".split(","):
                    df[name] = "Unclassified"

            df = df.groupby("kingdom")['percentage'].sum()
            # if a taxon is obsolete, the kingdom is empty. 
            # We will set the kingdom as Unclassified and raise a warning
            # if the count is > 5%
            if " " in df.index:
                percent = df.loc[" "]
                if percent > limit:
                    logger.warning("Found {}% of taxons with no kingdom".format(percent))
                if "Unclassified" in df.index:
                    df.loc['Unclassified'] += df.loc[' ']
                    df.drop(" ", inplace=True)
                else:
                    df.loc['Unclassified'] = df.loc[' ']
                    df.drop(" ", inplace=True)
            data[sample] = df

        df = pd.DataFrame(data)
        df = df.fillna(0)
        df = df.sort_index(ascending=False)
        df = df.sort_index(ascending=False, axis=1)

        #df.sum(axis=1).sort_values(ascending=False).index[0]

        return df

    def plot_stacked_hist(self, output_filename=None, dpi=200, kind="barh", 
        fontsize=10, edgecolor="k", lw=1, width=1, ytick_fontsize=10,
        max_labels=50):

        """Summary plot of reads classified.

        """
        df = self.get_df()

        fig, ax = pylab.subplots(figsize=(9.5,7))

        labels = []
        for kingdom in sorted(df.index):
            if kingdom == "Eukaryota":
                color="purple"
            elif kingdom == "Unclassified":
                color="grey"
            elif kingdom == "Bacteria":
                color="red"
            elif kingdom == "Viruses":
                color="darkgreen"
            elif kingdom == "Archaea":
                color="yellow"
            else:
                color="blue"
            if kind=="barh":
                pylab.barh(range(0, len(df.columns)), df.loc[kingdom], height=width,
                        left=df.loc[labels].sum().values, edgecolor=edgecolor,
                        lw=lw, color=color, alpha=0.8)
            else:
                pylab.bar(range(0, len(df.columns)), df.loc[kingdom], width=width,
                        bottom=df.loc[labels].sum().values, edgecolor=edgecolor,
                        lw=lw,color=color, alpha=0.8)
            labels.append(kingdom)


        positions = pylab.yticks()
        #ax.set_yticklabel(positions, labels, fontsize=ytick_fontsize)

        if kind == "barh":
            pylab.xlabel("Percentage (%)", fontsize=fontsize)
            pylab.ylabel("Sample index/name", fontsize=fontsize)
            if len(self.sample_names)<max_labels:
                pylab.yticks(range(len(self.sample_names)), self.sample_names[::-1], 
                    fontsize=ytick_fontsize)
            else:
                pylab.yticks([1], [""])
            pylab.xlim([0, 100])
            pylab.ylim([-0.5, len(df.columns)-0.5])
        else:
            pylab.ylabel("Percentage (%)", fontsize=fontsize)
            pylab.xlabel("Sample index/name", fontsize=fontsize)
            pylab.ylim([0, 100])
            if len(self.sample_names)<max_labels:
                pylab.xticks(range(len(self.sample_names)), self.sample_names[::-1], 
                    fontsize=ytick_fontsize)
            else:
                pylab.xticks([1], [""])
            pylab.xlim([-0.5, len(df.columns)-0.5])


        ax.legend(labels, title="kingdom",  bbox_to_anchor=(1,1))
        try:
            pylab.tight_layout()
        except:pass
        if output_filename:
            pylab.savefig(output_filename, dpi=dpi)




class MultiKrakenResults2():
    """Select several kraken output and creates summary plots

    ::

        import glob
        mkr = MultiKrakenResults2(glob.glob("*/*/summary.json"))
        mkr.plot_stacked_hist()

    """
    def __init__(self, filenames, sample_names=None):

        self.filenames = filenames
        if sample_names is None:
            self.sample_names = list(range(1,len(filenames)+1))
        else:
            self.sample_names = sample_names

    def get_df(self, limit=5, sorting_method="sample_name"):
        import pandas as pd
        import json
        data = {}
        for sample, filename in zip(self.sample_names, self.filenames):
            summary = json.loads(open(filename, "r").read())
            total = summary["total"]
            data[sample] = {
                    "unclassified": round(summary['unclassified']/total*100,2), 
                    "nreads": summary['total']}
            for db in summary['databases']:
                data[sample][db] = round(summary[db]['C'] / total * 100,2)
        df = pd.DataFrame(data)
        df = df.fillna(0)
        df = df.sort_index(ascending=False)
        df = df.sort_index(ascending=False, axis=1)


        # We sort the databases from best to worst. We ignore the nreads columns of
        # course, which is not related to the ability of each DB to classify the
        # reads
        index = df.loc[df.index.drop("nreads")].mean(axis=1).sort_values(ascending=False).index
        df = df.loc[list(index) + ["nreads"]]

        # sort the columns by sample names
        #df = df.sort_index(ascending=False, axis=1)
        #df.sum(axis=1).sort_values(ascending=False).index[0]
        return df

    def plot_stacked_hist(self, output_filename=None, dpi=200, kind="barh", 
            fontsize=10, edgecolor="k", lw=2, width=1, ytick_fontsize=10,
            max_labels=50, logx=False, alpha=0.8, colors=None, 
            cmap="hot_r", sorting_method="sample_name"):
        """Summary plot of reads classified.


        :param sorting_method: only by sample name for now
        """
        df = self.get_df(sorting_method=sorting_method)


        df = df.loc[["unclassified"]+[x for x in df.index if x!="unclassified"]]
        df = df.T
        del df['nreads']

        fig, ax = pylab.subplots(figsize=(9.5,7))

        labels = []
        # we will store grey/unclassified first and then other DB with a max of
        # 10 DBs
        L = len(df.columns) - 1
        from matplotlib import cm
        if colors is None:
            colors = [cm.get_cmap(cmap)(x) for x in pylab.linspace(0.2,1,L)]
        colors = ['grey'] + colors
        df.plot(kind="barh", stacked=True, width=width, 
                    edgecolor=edgecolor, color=colors, 
                    lw=lw,  alpha=alpha, ax=ax, legend=False)
        if logx is True:
            pylab.semilogx()


        positions = pylab.yticks()

        pylab.xlabel("Percentage (%)", fontsize=fontsize)
        pylab.ylabel("Sample index/name", fontsize=fontsize)
        if len(self.sample_names)<max_labels:
            pylab.yticks(range(len(self.sample_names)), df.index, 
                fontsize=ytick_fontsize)
        else:
            pylab.yticks([1], [""])
        pylab.xlim([0, 100])
        pylab.ylim([0-0.5, len(df.index)-0.5])  # +1 for the legends
        ax.legend(df.columns, title="DBs ",  bbox_to_anchor=(1,1))
        try:
            pylab.tight_layout()
        except:pass
        if output_filename:
            pylab.savefig(output_filename, dpi=dpi)
        return df
