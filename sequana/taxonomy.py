# -*- coding: utf-8 -*-
#
#  This file is part of Sequana software
#
#  Copyright (c) 2016-2020 - Sequana Development Team
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
import os
import re
from collections import Counter
from functools import wraps

from sequana import sequana_config_path
from sequana.lazy import pandas as pd
from sequana.utils.singleton import Singleton
from sequana.misc import wget

from easydev import TempFile

import colorlog
logger = colorlog.getLogger(__name__)


__all__ = ['NCBITaxonomy', 'Taxonomy']


class NCBITaxonomy():
    """


    """
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
            with TempFile() as fout_nodes:
                logger.info("Loading nodes.dmp from an URL {}".format(nodes))
                wget(nodes, fout_nodes.name)
                self.df_nodes = pd.read_csv(fout_nodes.name, sep="|", header=None)
        for i, _type in enumerate(self.df_nodes.dtypes):
            if _type == "O":
                self.df_nodes[i] = self.df_nodes[i].str.strip('\t')
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
            self.df_nodes.columns = ["taxid", "parent", "rank", 4, 5, "gc_id" ,"mt_id", 7, 8, 9, 10, 11, 12, 13]
            del self.df_nodes[13]
        except:
            self.df_nodes.columns = ["taxid", "parent", "rank", 4, 5]
            del self.df_nodes[5]


        # make sure they are ordered by taxon ID
        self.df_nodes.sort_values("taxid", inplace=True)
        self.df_nodes.set_index("taxid", inplace=True)

        # now we read the names
        if os.path.exists(names):
            self.df_names = pd.read_csv(names, sep="|", header=None)
        else:
            with TempFile() as fout_names:
                logger.info("Loading names.dmp from an URL {}".format(names))
                wget(names, fout_names.name)
                self.df_names = pd.read_csv(fout_names.name, sep="|", header=None)

        for i, _type in enumerate(self.df_names.dtypes):
            if _type == "O":
                self.df_names[i] = self.df_names[i].str.strip('\t')
        del self.df_names[4]
        self.df_names.columns = ['taxid', 'name', 'unique_name', 'key']
        self.df_names.set_index("taxid", inplace=True)

    def create_taxonomy_file(self, filename="taxonomy.dat"):
        logger.info("Please wait while creating the output file. "
            "This may take a few minutes")
        from easydev import Progress
        pb = Progress(len(self.df_nodes))
        count = 0
        df_names = self.df_names.query("key == 'scientific name'").copy()
        with open(filename, "w") as fout:

            for taxid in self.df_nodes.index:
                row = self.df_nodes.loc[taxid]
                fout.write("ID                        : {}\n".format(taxid))
                fout.write("PARENT ID                 : {}\n".format(row.parent))
                fout.write("RANK                      : {}\n".format(row['rank']))

                #names = df_names.loc[taxid]
                #print(
                fout.write("{:26s}: {}\n".format("SCIENTIFIC NAME", df_names.loc[taxid, "name"]))
                """    len(names)
                    for k,v in zip(names['key'], names['name']):
                        if k.upper() in ['SCIENTIFIC NAME', 'SYNONYM']:
                            fout.write("{:26s}: {}\n".format(k.upper(), v))
                except:
                    k, v = names['key'], names['name']
                    fout.write("{:26s}: {}\n".format(k.upper(), v))
                """
                fout.write("//\n")
                count += 1
                pb.animate(count)



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

    """
    def __init__(self, filename=None, verbose=True, online=True, source="ncbi"):
        """.. rubric:: constructor

        :param offline: if you do not have internet, the connction to Ensembl
            may hang for a while and fail. If so, set **offline** to True
        :param from: download taxonomy databases from ncbi
        """

        assert source in ['ncbi', 'ena']
        self.source = source

        if online:
            from bioservices import Ensembl, EUtils
            self.ensembl = Ensembl(verbose=False)

        self.records = {} # empty to start with.
        self.verbose = verbose

        if filename is None:
            self._dbname =  "taxonomy.dat"
            self.database = sequana_config_path + os.sep + self._dbname
        else:
            self.database = filename

        self._custom_db = sequana_config_path
        self._custom_db += "/taxonomy/taxonomy_custom.dat"

    def _update_custom_taxonomy_bases(self, taxid): # pragma: no cover
        """
        """
        taxid = str(taxid)
        self.eutils = EUtils(verbose=False)
        res = self.eutils.taxonomy_summary(taxid)
        if "error" in res[taxid]:
            print("not found in NCBI (EUtils)")
        else:
            print("found in NCBI (EUtils) and added to local databases")
            with open(self.custom_db, "w") as fout:
                data = res[taxid]
                fout.write("ID : {}\n".format(taxid))
                #fout.write("PARENT ID : {}\n".format(taxid))
                fout.write("RANK : {}\n".format(data['rank']))
                #fout.write("GC ID : {}\n".format(data['']))
                fout.write("SCIENTIFIC NAME : {}\n".format(data['scientificname']))

    def download_taxonomic_file(self, overwrite=False): #pragma: no cover
        """Loads entire flat file from EBI

        Do not overwrite the file by default.
        """
        import ftplib
        from sequana import sequana_config_path
        if os.path.exists(self.database) and overwrite is False:
            logger.info("Found taxonomy.dat file in sequana your path {}".format(sequana_config_path))
            return
        else:
            logger.info("Downloading and extracting the taxonomy file from the web. Please be patient.")

        if self.source == "ena":
            url = 'ftp.ebi.ac.uk'
        else:
            url = 'ftp.ncbi.nlm.nih.gov'

        self.ftp = ftplib.FTP(url)
        self.ftp.login()
        if self.source == "ena":
            # for the EBI ftp only: self.ftp.cwd('databases')
            self.ftp.cwd('pub')
            self.ftp.cwd('databases')
            self.ftp.cwd('taxonomy')
            logger.warning('Downloading and saving in %s. This is from ebi and may be behind the NCBI taxonomy' % self.database)
            self.ftp.retrbinary('RETR taxonomy.dat',
                open(self.database, 'wb').write)
            ftp.close()
        else:
            self.ftp.cwd('pub')
            self.ftp.cwd('taxonomy')
            logger.warning('Downloading and saving in %s from ncbi ftp' % self.database)
            import tempfile
            import shutil
            with tempfile.TemporaryDirectory() as tmpdir:
                filename = tmpdir + os.sep + "taxdump.tar.gz"
                self.ftp.retrbinary('RETR taxdump.tar.gz', open(filename, "wb").write)
                import tarfile
                tf = tarfile.open(filename)
                assert "nodes.dmp" in tf.getnames()
                assert "names.dmp" in tf.getnames()
                tf.extract("nodes.dmp", tmpdir)
                tf.extract("names.dmp", tmpdir)
                ncbi = NCBITaxonomy(tmpdir + os.sep + "names.dmp", tmpdir + os.sep + "nodes.dmp")
                ncbi.create_taxonomy_file(tmpdir  + os.sep + "taxonomy.dat")
                shutil.move(tmpdir + os.sep + "taxonomy.dat",
                            self.database)
            self.ftp.close()

    def load_records(self, overwrite=False):
        """Load a flat file and store records in :attr:`records`

        Since version 0.8.3 we use NCBI that is updated more often than the ebi
        ftp according to their README.

        ftp://ncbi.nlm.nih.gov/pub/taxonomy/

        """
        self.download_taxonomic_file(overwrite=overwrite)
        self.records = {}

        # TODO: check if it exists otherwise, load it ?
        if os.path.exists(self.database) is False:
            self.load()

        with open(self.database) as f:
            data = f.read().strip()

        # This is fast. tried parse package, much slower. cost of progress bar
        # is not important. 
        data = data.split("//\n") # the sep is //\n
        self._child_match = re.compile(r'ID\s+\:\s*(\d+)\s*')
        self._parent_match = re.compile(r'PARENT ID\s+\:\s*(\d+)\s*')
        self._rank_match = re.compile(r'RANK\s+\:\s*([^\n]+)\s*')
        self._name_match = re.compile(r'SCIENTIFIC NAME\s+\:\s*([^\n]+)\s*')

        from easydev import Progress
        pb = Progress(len(data))

        logger.info('Loading all taxon records.')
        for i, record in enumerate(data[0:]):
            dd = {'raw': record}
            dd['id'] = int(self._child_match.search(record).group(1))
            dd['parent'] = int(self._parent_match.search(record).group(1))
            dd['scientific_name'] = self._name_match.search(record).group(1)
            dd['rank'] = self._rank_match.search(record).group(1)
            self.records[dd["id"]] = dd
            if self.verbose:
                pb.animate(i+1)
        if self.verbose:
            print()
           
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
        taxon = int(taxon)
        lineage = self._gen_lineage_and_rank(taxon, [])
        lineage = [x[0] for x in lineage]
        return lineage

    @load_taxons
    def _gen_lineage_and_rank(self, taxon, lineage_rank=[]):
        # recursively filling the lineage argument

        try:
            record = self.records[taxon]
        except:
            return [('unknown_taxon:{}'.format(taxon), 'no rank')]

        parent = int(record['parent'])

        if taxon == 1:
            lineage_rank.append((record['scientific_name'], record['rank']))
            lineage_rank.reverse()
            return lineage_rank
        else:
            lineage_rank.append((record['scientific_name'], record['rank']))
            return self._gen_lineage_and_rank(parent, lineage_rank)

    @load_taxons
    def get_parent_taxon(self, taxon):
        return self.records[taxon]['parent']

    @load_taxons
    def get_parent_name(self, taxon):
        taxid = self.get_parent_taxon(taxon)
        return self.records[taxid]['scientific_name']

    @load_taxons
    def get_lineage_and_rank(self, taxon):
        """Get lineage and rank of a taxon

        :param int taxon:
        :return: a list of tuples. Each tuple is a pair of taxon name/rank
            The list is the lineage for to the input taxon.

        """
        taxon = int(taxon)
        lineage = self._gen_lineage_and_rank(taxon, [])
        return lineage

    @load_taxons
    def get_ranks(self):
        return  Counter([x['rank'] for x in self.records.values()])

    @load_taxons
    def get_record_for_given_rank(self, rank):
        return [x for x in self.records.values() if x['rank'] == rank]

    @load_taxons
    def get_names_for_given_rank(self, rank):
        data = [x for x in self.records.values() if x['rank'] == rank]
        return [x['scientific_name'] for x in data]

    @load_taxons
    def get_children(self, taxon):
        taxon = str(taxon)
        children = [self.records[k] for k in self.records.keys()
                if self.records[k]['parent'] == taxon]
        children = [child['id'] for child in children]
        return children

    @load_taxons
    def __getitem__(self, iden):
        return self.records[iden]

    @load_taxons
    def __len__(self):
        return len(self.records)

    def append_existing_database(self, filename):
        """

        Taxonomy DB looks like::

            ID                        : 2731450
            PARENT ID                 : 1914233
            RANK                      : genus
            SCIENTIFIC NAME           : Limnoglobus
            //


            a = NCBITaxonomy("names.dmp", "nodes.dmp")
            a.create_taxonomy_file("taxonomy.dat")
            tax = Taxonomy()
            tax.append_existing_database("taxonomy.dat")
        """
        tax = Taxonomy(filename)
        tax.load_records()
        self.load_records()
        toadd = []
        for record in tax.records.keys():
            if record not in self.records:
                toadd.append(record)

        with open(self.database, "a") as fout:
            for record in toadd:
                fout.write(tax.records[record]['raw']+"//\n")
