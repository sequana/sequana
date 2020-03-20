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
import collections
import os
import ftplib
import subprocess
import sys
import glob


from sequana.databases import ENADownload
from easydev import execute, TempFile, Progress, md5, DevTools

from sequana.lazy import pandas as pd
from sequana import logger
logger.name = __name__



__all__ = ["KrakenBuilder"]


class KrakenBuilder():
    """This class will help you building a custom Kraken database


    You will need a few steps, and depending on the FASTA files you want to
    include lots of resources (memory and space wise). In the following example,
    we will be reasonable and use only viruses FASTA files.

    First, we need to create the data structure directory. Let us call it
    **virusdb**::

        from sequana import KrakenBuilder
        kb = KrakenBuilder("virusdb")

    We then need to download a large taxonomic database from NCBI. You may
    already have a local copy, in which case you would need to copy it in
    virusdb/taxonomy directory. If not, type::

        kb.download_taxonomy()

    The virusdb/taxonomy directory will contain about 8.5G of data.

    Note that this currently requires the unix tools **wget** and **tar**.

    Then, we need to add some fasta files. You may download specific FASTA files
    if you know the accession numbers using :meth:`download_accession`. However,
    we also provide a method to download all viruses from ENA::

        kb.download_viruses()

    This will take a while to download the more than 4500 FASTA files (10
    minutes on a good connection). You will end up with a data set of about 100
    Mb of FASTA files.

    If you wish to download other FASTA (e.g. all bacteria), you will need to
    use another class from the :mod:`sequana.databases`::

        from sequana.databases import ENADownload
        ena = ENADownload()
        ena.download_fasta("bacteria.txt", output_dir="virusdb/library/added")

    Please see the documentation for more options and list of species to
    download.

    It is now time to build the DB itself. This is based on the kraken tool.
    You may do it yourself in a shell::

        kraken-build  --rebuild -db virusdb --minimizer-len 10 --max-db-size 4 --threads 4
        --kmer-len 26 --jellyfish-hash-size 500000000

    Or you the KrakenBuilder. First you need to look at the :attr:`params`
    attribute. The most important key/value that affect the size of the DB are::

        kb.params['kmer_length']  (max value is 31)
        kb.params['max_db_size'] is tha max size of the DB files in Gb
        kb.params['minimizer_len']

    To create a small DB quickly, we set those values::

        kb.params['kmer_length']  = 26
        kb.params['minimizer_len'] = 10

    However, for production, we would recommend 31 and 13 (default)

    This takes about 2 minutes to build and the final DB is about 800Mb.

    Lots of useless files are in the direcory and can be removed using kraken
    itself. However we do a little bit more and therefore have our own
    cleaning function::

        kb.clean_db()

    Kraken-build uses jellyfish. The **hash_size** parameter is the jellyfish
    hash_size parameter. If you set it to 6400M, the memory required is about
    6.9bytes times 6400M that is 40Gb of memory. The default value used here
    means 3.5Gb are required.

    The size to store the DB itself should be

    :math:

        sD + 8 (4^M)

    where **s** is about 12 bytes (used to store a kmer/taxon pair, D is the
    number of kmer in the final database, which cannot be estimated before
    hand, and M the length minimiser parameter.


    The quick way:
    =====================

        kb = KrakenBuilder("virusdb")
        kb.run(['virus']) # use only viruses from ENA list

    Here, you may want to re-run the analysis with different parameters
    for the database built. If you require the virus DB, it has been
    downloaded already so this step will be skip. The Taxon DB does not
    need to be downloaded again, so set download_taxonomy to False.

    Before, let us change the parameter to build a full database::

        kb.params['kmer_length']  = 31
        kb.params['minimizer_len'] = 13

    We have here instead of 800Mb DB a new DB of 1.5Gb but it should
    take more or less the same time to build it

    Finally if you do not need to test it anymore, you may clean the DB once for
    all. This will remove useless files. The directory's name is the name of the
    DB that should be used in e.g. the quality_control pipeline. To clean the
    data directory, type::

        kb.clean_db()

    """
    def __init__(self, dbname):
        """.. rubric:: Constructor

        :param str dbname: Create the Kraken DB in this directory


        """
        # See databases.py module
        self.dbname = dbname
        self.enadb = ENADownload()
        self.valid_dbs = self.enadb._metadata.keys()

        # mini_kraken uses minimiser-length = 13, max_db =4, others=default so
        # kmer-len=31 hashsize=default
        self.params = {
            "dbname": self.dbname,
            "minimizer_len": 10,
            "max_db_size": 4,
            "threads": 4,
            "kmer_length" : 26,
            "hash_size" : 500000000
        }

        self.init()

    def init(self):
        # mkdir library
        self.library_path = self.dbname + os.sep + "library"
        self.taxon_path = self.dbname + os.sep + "taxonomy"
        self.fasta_path = self.library_path + os.sep + "added"

        self._devtools = DevTools()
        self._devtools.mkdir(self.dbname)
        self._devtools.mkdir(self.library_path)
        self._devtools.mkdir(self.fasta_path)
        self._devtools.mkdir(self.taxon_path)

    def download_ncbi_refseq(self, category):
        """Download all files of type *fna* from ncbi FTP.

        ::

            kb = KrakenBuilder()
            kb.download_ncbi_refseq("viral")

        """
        import ftplib
        import os

        ftp = ftplib.FTP("ftp.ncbi.nlm.nih.gov")
        ftp.login("anonymous", "thomas.cokelaer@pasteur.fr")
        ftp.cwd("refseq/release/{}".format(category))

        import io
        file_mapper = {}
        for filename in ftp.nlst():
            if "genomic.fna" in filename:
                ftp.retrbinary('RETR ' + filename , open(filename, "wb").write)
                print(filename)

    def download_accession(self, acc):
        """Donwload a specific Fasta from ENA given its accession number

        Note that if you want to add specific FASTA from ENA, you must use
        that function to make sure the header will be understood by Kraken;
        The header must use a GI number (not ENA)

        """
        output = self.dbname+os.sep + "library" + os.sep + "added"
        """Download a specific FASTA file given its ENA accession number """
        self.enadb.download_accession(acc, output=output)

    def download_viruses(self):
        self.enadb.download_fasta("virus.txt", output_dir=self.fasta_path)

    def run(self, dbs=[], download_taxon=True):
        """Create the Custom Kraken DB

        #. download taxonomy files
        #. Load the DBs (e.g. viruses)
        #. Build DB with kraken-build
        #. Clean it up

        """
        # Start with the FASTA
        self._download_dbs(dbs)

        self.download_taxonomy()

        # search for taxon file. If not found, error
        required = self.taxon_path + os.sep + "gi_taxid_nucl.dmp"

        if required not  in glob.glob(self.taxon_path + os.sep + "*"):
            raise IOError("Taxon file not found")

        print("\nDepending on the input, this step may take a few hours to finish")
        self._build_kraken()

    def download_taxonomy(self, force=False):
        """Download kraken data

        The downloaded file is large (1.3Gb) and the unzipped file is about 9Gb.

        If already present, do not download the file except if the *force*
        parameter is set to True.

        """

        # If the requested file exists, nothing to do
        expected_filename = self.taxon_path + os.sep + "gi_taxid_nucl.dmp"
        expected_md5 = "8c182ac2df452d836206ad13275cd8af"
        print('\nDownloading taxonomy files. Takes a while depending on your connection')

        if os.path.exists(expected_filename) is False or \
                md5(expected_filename) != expected_md5:
            # download taxonomy
            # We could use kraken-build --download-taxonomy + a subprocess but
            # even simpler to get the file via ftp
            FTP = "ftp.ncbi.nih.gov"
            execute("wget %s/pub/taxonomy/gi_taxid_nucl.dmp.gz --directory-prefix %s"
                % (FTP, self.taxon_path))
            # Unzip the files
            execute('unpigz %s/gi_taxid_nucl.dmp.gz' % self.taxon_path)
        else:
            print("Found local expected file %s " % expected_filename)

        expected_filename = self.taxon_path + os.sep + "names.dmp"
        expected_md5 = "90d88912ad4c94f6ac07dfab0443da9b"
        if os.path.exists(expected_filename) is False or \
                md5(expected_filename) != expected_md5:

            execute("wget %s/pub/taxonomy/taxdump.tar.gz --directory-prefix %s"
                % (FTP, self.taxon_path))

            execute('tar xvfz %s/taxdump.tar.gz -C %s' %
                (self.taxon_path, self.taxon_path))
        else:
            print("Found local expected file %s " % expected_filename)

    def _download_dbs(self, dbs=[]):
        print("Downloading all Fasta files for %s" % dbs)
        # Download the DBs in it
        from .databases import ENADownload
        for db in dbs:
            if db not in self.valid_dbs and os.path.exists(db) is False:
                msg = "db must be a local file with a list of ENA or one of"
                for this in self.ena._metadata.keys():
                    msg += " - %s" % this
                raise ValueError(msg)
            self.ena.download_fasta(db, output_dir=self.fasta_path)

    def _build_kraken(self):
        print('Building the kraken db ')
        self.params['hash_size'] = int(self.params["hash_size"])

        cmd = """kraken-build  --rebuild -db %(dbname)s \
            --minimizer-len %(minimizer_len)s\
            --max-db-size %(max_db_size)s \
            --threads %(threads)s\
            --kmer-len %(kmer_length)s \
            --jellyfish-hash-size %(hash_size)s""" % self.params

        # again, kraken-build prints on stderr so we cannot use easydev.shellcmd
        execute(cmd)

    def clean_db(self):
        """Once called, you will not be able to append more FASTA files

        """
        # Now we can clean the kraken db:
        print('Cleaning the kraken db ')
        # Clean the nodes.dmp and names.dmp
        print('Identifying the GI numbers')
        gis = self.get_gis()
        taxons = self.get_taxons_from_gis(gis)
        print("")

        self.gis = gis
        self.taxons = taxons

        # This cleans the nodes.dmp and names.dmp. This must be done
        # before kraken-build --clean since it requires the gi_taxid_nucl.dmp
        # file
        names_file = self.taxon_path + os.sep + "names.dmp"
        nodes_file = self.taxon_path + os.sep + "nodes.dmp"
        names_file_temp = self.taxon_path + os.sep + "names_temp.dmp"
        nodes_file_temp = self.taxon_path + os.sep + "nodes_temp.dmp"

        taxon_file_reader = NCBITaxonReader(names=names_file, nodes=nodes_file,
            verbose=True)
        print("Filtering")
        taxon_file_reader.filter_nodes_dmp_file(nodes_file, nodes_file_temp,
            taxons=taxons)
        taxon_file_reader.filter_names_dmp_file(names_file, names_file_temp,
            taxons=taxons)

        # mv the new files into the old ones
        os.rename(names_file_temp, names_file)
        os.rename(nodes_file_temp, nodes_file)

        # Finally, the kraken cleaning itself
        cmd = "kraken-build --clean --db %s" % self.params['dbname']
        execute(cmd)

    def get_gis(self, extensions=['fa']):
        self.filenames = []
        root = self.dbname
        for extension in extensions:
            self.filenames.extend( list(glob.iglob("%s/library/**/*%s" %
                (root, extension))))
        for extension in extensions:
            self.filenames.extend( list(glob.iglob("%s/library/**/**/*%s" %
                (root, extension))))

        N = len(self.filenames)
        pb = Progress(N)
        gis = []
        for i, filename in enumerate(self.filenames):
            data = open(filename, "r")
            line = data.readline()
            if line.startswith('>'):
                assert "gi" in line, "expected >gi to be found at the beginning"
                gi = line[1:].split("|")[1]
            else:
                raise ValueError("This file %s does not seem to be a FASTA file" % filename)
            gis.append(gi)
            pb.animate(i+1)
        print()
        gis = [int(x) for x in gis]
        self.gis = gis

        assert len(gis) == len(self.filenames)
        return gis

    def get_taxons_from_gis(self, gis, filename="gi_taxid_nucl.dmp"):
        filename = self.taxon_path + os.sep + filename
        data = pd.read_csv(filename, chunksize=1000000, sep='\t',
            header=None)
        N = 560  # with time this number will be deprecated but good for now

        local_gis = gis[:]

        # We will found GI an order than different from the input gis list so
        # we will need to keep track of the order
        found_gis = []
        taxons = [32644] * len(gis)   # 32644 means unidentified
        # we search for the unique gis. Once found, we remove them from the
        # vector and keep going until the vector is empty or there is no more
        # chunks. A good sanity check is that the final gis vector should be
        # empty meaning all have been found. We do not care about the order
        # of the final taxons vector as compare to the GI vector

        print("Scanning %s to look for %s GI numbers" % (filename, len(gis)))
        pb = Progress(N)
        for i, chunk in enumerate(data):
            chunk.set_index(0, inplace=True)
            chunk = chunk.loc[local_gis].dropna()

            # keep the GI and Taxon
            found_gis.extend([int(x) for x in list(chunk.index)])
            
            # update the remaining GIs and the taxons
            for gi, tax in zip(chunk.index, chunk.values):
                local_gis.remove(gi)
                index = gis.index(gi)
                taxons[index] = tax

            # no need to carry on if all GIs were found
            if len(local_gis) == 0:
                break
            pb.animate(i+1)
        print("")

        taxons = [int(x) for x in taxons]
        return taxons


class NCBITaxonReader(object):
    """This class will help in reading, handling, simplifying NCBI taxon DB
    used by Kraken

    When downloading NCBI taxonomy DB using e.g. Kraken, we end up with very
    large files. One is called names.dmp and the other nodes.dmp.

    Yet, when we build a Kraken DB, we select a subset of fasta/species. It
    would be convenient to be able to filter the names.dmp and nodes.dmp to reduce
    the size of the final Kraken DB.

    The names.dmp is just a CSV file. The header looks like::


        1   |   all |       |   synonym |
        1   |   root    |       |   scientific name |
        2   |   Bacteria    |   Bacteria <prokaryote>   |   scientific name |
        2   |   Monera  |   Monera <Bacteria>   |   in-part |
        2   |   Procaryotae |   Procaryotae <Bacteria>  |   in-part |

    It is a tabulated file. If we ignore the | signs, it contains 4 columns::

        taxid
        name
        unique name
        type of name

    The *unique name* column is generally empty and is dropped internally.
    There are different types of *name*, so there can be several rows for
    a given *taxid*. For instance
    for the taxon 1, there isa  *scientific name* and a **synonym**

    The :attr:`df_name` is a dataframe that stores the taxid, name and type of
    name in a dataframe.

    The second file 'nodes.dmp') looks like::

        1 | 1       | no rank |     | 8 | 0 | 1  | 0  | 0 | 0 | 0 | 0 |   |
        2 | 131567  | superkingdom  |   | 0 | 0  | 11 | 0 | 0 | 0 | 0 | 0 | |
        6 | 335928  | genus   |     | 0 | 1 | 11 | 1  | 0 | 1 | 0 | 0 |   |
        7 | 6       | species | AC  | 0 | 1 | 11 | 1  | 0 | 1 | 1 | 0 |   |
        9 | 32199   | species | BA  | 0 | 1 | 11 | 1  | 0 | 1 | 1 | 0 |   |

    Again this is a tabulated file. The first three columns are taxid, parent taxid,
    and rank. Rank is species, genus, family, phylum, etc.

    ::

        gi = kraken.GetGIFromLibrary("../library/")
        gis = gi.get_gis()
        k = kraken.CreateGI_to_TaxonFile()
        taxons = k.get_taxons_from_gis(gis)

        kk = kraken.NCBITaxonReader()
        kk.filter_nodes_dmp_file(taxons=taxons)
        kk.filter_names_dmp_file(taxons=taxons)

    """
    def __init__(self, names="names.dmp", nodes="nodes.dmp", verbose=True):
        if verbose:
            print("Reading %s" %  names)
        self.df_name = pd.read_csv(names, sep='\t', header=None)
        self.df_name = self.df_name[[0,2,6]]
        self.df_name.columns = ["taxon", "name", "scname"]

        # This will provide a faster lookup table to search for scientific
        # names given a taxon. We can drop rows that are not scientific names
        # and set the taxons as index
        _subdf = self.df_name.query("'scientific name' in scname")
        self._subdf = _subdf.set_index("taxon")

        # Here, this is for general purpose (slower it we were to use
        # this for the get_scientic_name method
        self._group_name = self.df_name.groupby('taxon').groups

        if verbose:
            print("Reading %s" %  nodes)
        self.df_nodes = pd.read_csv(nodes, sep='\t', header=None)
        self.df_nodes = self.df_nodes[[0,2,4,6,8,10,12,14,16,18,20,22,24]]
        self.df_nodes.columns = ['taxon', 'parent_taxon', 'rank', 0,1,2,3,4,5,6,7,8,9]

        self._df_nodes_taxon = self.df_nodes.copy()
        self._df_nodes_taxon.set_index('taxon', inplace=True)

    def get_number_taxon(self):
        """Return number of unique taxon"""
        return len(self.df_name['taxon'].unique())

    def get_average_name_per_taxon(self):
        """Return number of rows/names per node/taxon"""
        return pylab.mean([len(values) for values in self._group_name.values()])

    def get_scientific_name(self, taxon):
        """Return scientific name of a given Taxon"""
        # Takes 2 minutes to scan all taxons
        return self._subdf.loc[taxon].values[0]

    def get_taxon_from_scientific_name(self, scname):
        """Return taxon corresponding to a scientific name

        return: unique taxon or first one found. If none found, returns None
        """
        res = self.df_name.query("@scname in name")['taxon']
        return res

    def search(self, name):
        """Search names colum"""
        return self.df_name[self.df_name['name'].apply(lambda x : name in x)]

    def get_family(self, taxon):
        """Get all parent taxons"""
        taxons = [1]
        df = self._df_nodes_taxon
        while True:
            res = df.loc[taxon]
            taxons.append(taxon)
            taxon = res['parent_taxon']
            # hopefully there is always a family link to 0 or 1
            if len(res) == 0 or taxon in [0,1]:
                break
        return taxons

    def filter_nodes_dmp_file(self, filename="nodes.dmp",
            output="nodes_filtered.dmp", taxons=[]):
        with open(filename, "r") as fin:
            with open(output, "w") as fout:
                for line in fin.readlines():
                    if int(line.split("\t", 1)[0]) in taxons:
                        fout.write(line)

    def filter_names_dmp_file(self, filename="names.dmp",
            output="names_filtered.dmp", taxons=[]):

        all_taxons = set()
        pb = Progress(len(taxons))
        for i, taxon in enumerate(taxons):
            parents = self.get_family(taxon)
            all_taxons.update(parents)
            pb.animate(i+1)
        print("")
        with open(filename, "r") as fin:
            with open(output, "w") as fout:
                for line in fin.readlines():
                    if int(line.split("\t", 1)[0]) in all_taxons:
                        fout.write(line)

