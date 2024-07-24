#
#  This file is part of Sequana software
#
#  Copyright (c) 2016-2020 - Sequana Development Team
#
#  Distributed under the terms of the 3-clause BSD license.
#  The full license is in the LICENSE file, distributed with this software.
#
#  website: https://github.com/sequana/sequana
#  documentation: http://sequana.readthedocs.io
#
##############################################################################
import ftplib
import gzip
import os
import shutil
import tarfile
import tempfile
from functools import wraps
from pathlib import Path

import colorlog
from tqdm import tqdm

from sequana import sequana_config_path
from sequana.lazy import pandas as pd
from sequana.misc import wget
from sequana.utils.singleton import Singleton

logger = colorlog.getLogger(__name__)


__all__ = ["NCBITaxonomy", "Taxonomy"]


class NCBITaxonomy:
    """ """

    def __init__(self, names, nodes):
        """

        :param names: can be a local file or URL
        :param nodes: can be a local file or URL

        """
        # Path to existing files
        logger.info("Reading input files")
        self.names = names
        self.nodes = nodes

        # First, the nodes
        if os.path.exists(nodes):
            self.df_nodes = pd.read_csv(nodes, sep="|", header=None)
        else:
            with tempfile.TemporaryFile() as fout_nodes:
                logger.info("Loading nodes.dmp from an URL {}".format(nodes))
                wget(nodes, fout_nodes.name)
                self.df_nodes = pd.read_csv(fout_nodes.name, sep="|", header=None)

        for i, _type in enumerate(self.df_nodes.dtypes):
            if _type == "O":
                self.df_nodes[i] = self.df_nodes[i].str.strip("\t")
        """
        tax_id                  -- node id in GenBank taxonomy database
        parent tax_id               -- parent node id in GenBank taxonomy database
        rank                    -- rank of this node (superkingdom, kingdom, ...)
        embl code               -- locus-name prefix; not unique
        division id             -- see division.dmp file
        inherited div flag  (1 or 0)        -- 1 if node inherits division from parent
        genetic code id             -- see gencode.dmp file
        inherited GC  flag  (1 or 0)        -- 1 if node inherits genetic code from parent
        mitochondrial genetic code id       -- see gencode.dmp file
        inherited MGC flag  (1 or 0)        -- 1 if node inherits mitochondrial gencode from parent
        GenBank hidden flag (1 or 0)            -- 1 if name is suppressed in GenBank entry lineage
        hidden subtree root flag (1 or 0)       -- 1 if this subtree has no sequence data yet
        comments                -- free-text comments and citations
        """

        try:
            self.df_nodes.columns = [
                "taxid",
                "parent",
                "rank",
                4,
                5,
                "gc_id",
                "mt_id",
                7,
                8,
                9,
                10,
                11,
                12,
                13,
            ]
            del self.df_nodes[13]
        except Exception:  # pragma: no cover
            self.df_nodes.columns = ["taxid", "parent", "rank", 4, 5]
            del self.df_nodes[5]

        # make sure they are ordered by taxon ID
        self.df_nodes.sort_values("taxid", inplace=True)
        self.df_nodes.set_index("taxid", inplace=True)

        # now we read the names
        if os.path.exists(names):
            self.df_names = pd.read_csv(names, sep="|", header=None)
        else:
            with tempfile.TemporaryFile() as fout_names:
                logger.info(f"Loading names.dmp from an URL {names}")
                wget(names, fout_names.name)
                self.df_names = pd.read_csv(fout_names.name, sep="|", header=None)

        for i, _type in enumerate(self.df_names.dtypes):
            if _type == "O":
                self.df_names[i] = self.df_names[i].str.strip("\t")
        del self.df_names[4]
        self.df_names.columns = ["taxid", "name", "unique_name", "key"]
        self.df_names.set_index("taxid", inplace=True)

    def create_taxonomy_file(self, filename="taxonomy.csv.gz"):

        filename = Path(filename)
        logger.info("Please wait while creating the output file. This may take a few minutes")
        if filename.suffixes != [".csv", ".gz"]:
            raise ValueError(f"{filename} extension must be '.csv.gz'")

        count = 0
        df_names = self.df_names.query("key == 'scientific name'").copy()

        # first we create the CSV file
        logger.info("Creating CSV file")
        with filename.with_suffix("").open(mode="w") as fout:
            fout.write("id,parent,rank,scientific_name\n")

            for taxid in tqdm(self.df_nodes.index):
                row = self.df_nodes.loc[taxid]
                sc = df_names.loc[taxid, "name"]
                fout.write(f"{taxid},{row.parent},{row['rank']},{sc}\n")

        # second, we gzip it. input must be read as binary
        logger.info("Compressing CSV file")
        with filename.with_suffix("").open(mode="rb") as fin:
            bindata = fin.read()
            with gzip.open(filename, "wb") as f:
                f.write(bindata)


def load_taxons(f):
    @wraps(f)
    def wrapper(*args, **kargs):
        if len(args[0].records) == 0:
            args[0].load_records()

        return f(*args, **kargs)

    return wrapper


class Taxonomy(metaclass=Singleton):
    """This class should ease the retrieval and manipulation of Taxons

    There are many resources to retrieve information about a Taxon.
    For instance, from BioServices, one can use UniProt, Ensembl, or
    EUtils. This is convenient to retrieve a Taxon (see :meth:`fetch_by_name`
    and :meth:`fetch_by_id` that rely on Ensembl). However, you can
    also download a flat file from EBI ftp server, which
    stores a set or records (2.8M (april 2020).

    Note that the Ensembl database does not seem to be as up to date
    as the flat files but entries contain more information.

    for instance taxon 2 is in the flat file but not available through
    the :meth:`fetch_by_id`, which uses ensembl.

    So, you may access to a taxon in 2 different ways getting differnt
    dictionary. However, 3 keys are common (id, parent, scientific_name)

    ::

        >>> t = taxonomy.Taxonomy()
        >>> t.fetch_by_id(9606)   # Get a dictionary from Ensembl
        >>> t.records[9606] # or just try with the get
        >>> t[9606]
        >>> t.get_lineage(9606)


    Possible ranks are various. You may have biotype, clade, etc ub generally speaking
    ranks are about lineage. For a given rank, e.g. kingdom, you may have sub division such
    as superkingdom and subkingdom. order has even more subdivisions (infra, parv, sub, super)

    Since version 0.8.3 we use NCBI that is updated more often than the ebi
    ftp according to their README. ftp://ncbi.nlm.nih.gov/pub/taxonomy/
    We use Ensembl to retrieve various information regarding taxons.
    """

    def __init__(self, filename=None, verbose=True, online=True, source="ncbi"):
        """.. rubric:: constructor

        :param offline: if you do not have internet, the connction to Ensembl
            may hang for a while and fail. If so, set **offline** to True
        :param from: download taxonomy databases from ncbi
        """
        assert source in ["ncbi", "ena"]
        self.source = source

        if online:
            from bioservices import Ensembl, EUtils

            self.ensembl = Ensembl(verbose=False)

        self.records = {}  # empty to start with.
        self.verbose = verbose

        if filename is None:
            self._dbname = "taxonomy.csv.gz"
            self.database = sequana_config_path + os.sep + self._dbname
        else:
            assert str(filename).endswith(".csv.gz")
            self.database = filename

        self._custom_db = sequana_config_path
        self._custom_db += "/taxonomy/taxonomy_custom.csv.gz"

    def download_taxonomic_file(self, overwrite=False):  # pragma: no cover
        """Loads entire flat file from NCBI

        Do not overwrite the file by default.
        """

        if os.path.exists(self.database) and overwrite is False:
            logger.info(f"Found {self.database} file in sequana your path {sequana_config_path}")
            return

        if self.source == "ena":
            url = "ftp.ebi.ac.uk"
        else:
            url = "ftp.ncbi.nlm.nih.gov"

        self.ftp = ftplib.FTP(url)
        self.ftp.login()
        if self.source == "ena":
            # for the EBI ftp only: self.ftp.cwd('databases')
            self.ftp.cwd("pub")
            self.ftp.cwd("databases")
            self.ftp.cwd("taxonomy")
            logger.warning("Downloading and saving data in {self.database}")
            self.ftp.retrbinary("RETR taxonomy.dat", open(self.database, "wb").write)
            self.ftp.close()
        else:
            self.ftp.cwd("pub")
            self.ftp.cwd("taxonomy")
            logger.warning(f"Downloading and saving data in {self.database}")

            with tempfile.TemporaryDirectory() as tmpdir:
                filename = tmpdir + os.sep + "taxdump.tar.gz"
                self.ftp.retrbinary("RETR taxdump.tar.gz", open(filename, "wb").write)

                tf = tarfile.open(filename)
                assert "nodes.dmp" in tf.getnames()
                assert "names.dmp" in tf.getnames()
                logger.info("Extracting nodes.dmp")
                tf.extract("nodes.dmp", tmpdir)
                logger.info("Extracting names.dmp")
                tf.extract("names.dmp", tmpdir)

                logger.info("Extracting and building local database")
                ncbi = NCBITaxonomy(tmpdir + os.sep + "names.dmp", tmpdir + os.sep + "nodes.dmp")

                tmp_taxfile = tmpdir + os.sep + "taxonomy.csv.gz"
                ncbi.create_taxonomy_file(tmp_taxfile)
                shutil.move(tmp_taxfile, self.database)

            self.ftp.close()

    def load_records(self, overwrite=False):
        """Load a flat file and store records in :attr:`records`"""
        self.download_taxonomic_file(overwrite=overwrite)

        # for usecols=range(4) cause last column may contain extra commas
        try:
            self.records = pd.read_csv(self.database, index_col=0, compression="gzip", usecols=range(4))
        except gzip.BadGzipFile:
            logger.error(f"input file {self.database} should be gzipped")
            raise gzip.BadGzipFile

    def load_records_from_csv(self, filename):
        df = pd.read_csv(filename)

    def find_taxon(self, taxid, mode="ncbi"):
        taxid = str(taxid)
        if mode == "ncbi":
            from bioservices import EUtils

            self.eutils = EUtils(verbose=False)
            res = self.eutils.taxonomy_summary(taxid)
        else:
            res = self.ensembl.get_taxonomy_by_id(taxid)
        return res

    @load_taxons
    def fetch_by_id(self, taxon):
        """Search for a taxon by identifier

        :return; a dictionary.

        ::

            >>> ret = s.search_by_id('10090')
            >>> ret['name']
            'Mus Musculus'

        """
        res = self.ensembl.get_taxonomy_by_id(taxon)
        return res

    @load_taxons
    def fetch_by_name(self, name):
        """Search a taxon by its name.

        :param str name: name of an organism. SQL cards possible e.g.,
            _ and % characters.
        :return: a list of possible matches. Each item being a dictionary.

        ::

            >>> ret = s.search_by_name('Mus Musculus')
            >>> ret[0]['id']
            10090

        """
        res = self.ensembl.get_taxonomy_by_name(name)
        return res

    @load_taxons
    def get_lineage(self, taxon):
        """Get lineage of a taxon

        :param int taxon: a known taxon
        :return: list containing the lineage

        """
        # important to reinit the second argument to []
        lineage = self._gen_lineage_and_rank(taxon, [])
        return [x[0] for x in lineage]

    @load_taxons
    def _gen_lineage_and_rank(self, taxon, lineage_rank=[]):
        # recursively filling the lineage argument

        try:
            record = self.records.loc[taxon]
        except Exception as err:
            print(err)
            return [("unknown_taxon:{}".format(taxon), "no rank")]
        parent = int(record["parent"])

        if taxon == 1:
            lineage_rank.append((record["scientific_name"], record["rank"]))
            lineage_rank.reverse()
            return lineage_rank
        else:
            lineage_rank.append((record["scientific_name"], record["rank"]))
            return self._gen_lineage_and_rank(parent, lineage_rank)

    @load_taxons
    def get_parent_taxon(self, taxon):
        return self.records.loc[taxon, "parent"]

    @load_taxons
    def get_parent_name(self, taxon):
        taxid = self.get_parent_taxon(taxon)
        return self.records.loc[taxid, "scientific_name"]

    @load_taxons
    def get_lineage_and_rank(self, taxon):
        """Get lineage and rank of a taxon

        :param int taxon:
        :return: a list of tuples. Each tuple is a pair of taxon name/rank
            The list is the lineage for to the input taxon.

        """
        return self._gen_lineage_and_rank(taxon, [])

    @load_taxons
    def get_ranks(self):
        return self.records.groupby("rank").count().parent.sort_values()

    def get_records_for_given_rank(self, rank):
        return self.records.query("rank == @rank")

    @load_taxons
    def get_names_for_given_rank(self, rank):
        return self.records.query("rank == @rank").scientific_name.unique()

    @load_taxons
    def get_children(self, taxon):
        return self.records.query("parent == @taxon").index.values

    @load_taxons
    def __len__(self):
        return len(self.records)

    def __getitem__(self, taxon):
        return self.records.loc[taxon]
