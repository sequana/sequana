#
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
"""Utilities to access to online FASTA, taxon, lineage ..."""
import os
import glob
import math
import ftplib

from easydev import AttrDict, execute, Progress

from sequana.lazy import pandas as pd

from bioservices import ENA, EUtils

import colorlog

logger = colorlog.getLogger(__name__)


class EUtilsTools(object):
    """Utilities to fetch information about accession numbers

    ::

        >>> from sequana.databases import EUtilsTools
        >>> et = EUtilsTools()
        >>> et.accession_to_info("K01711.1")
        {'K01711.1': {'accession': '331784',
          'comment': 'Measles virus (strain Edmonston), complete genome',
          'gi': '331784',
          'identifier': 'gi|331784|gb|K01711.1|MEANPCG[331784]',
          'taxid': '11234'}}

    """

    def __init__(self):
        self.eutils = EUtils(cache=True)

    def get_fasta(self, accession):
        data = self.eutils.EFetch(
            "nucleotide", id=accession, retmode="text", rettype="fasta"
        )
        return data

    def accession_to_info(self, ids):
        """An accession or list of them returns list of dictionaries"""
        res = self.eutils.EFetch(db="nuccore", id=ids, rettype="docsum", retmode="json")
        res = res["result"]

        # now we can loop over all identifiers
        records = {}

        for i, uid in enumerate(res["uids"]):
            data = res[uid]
            accession = data["accessionversion"]

            record = {
                "taxid": data["taxid"],
                "accession": data["accessionversion"],
                "gi": data["gi"],
                "comment": data["title"],
            }

            records[accession] = AttrDict(**record)
        return records


class ENADownload(object):
    """Downloader to retrieve genome fasta files from ENA amongst other things

    In order to facilitate the download of FASTA files (e.g. to build a Kraken
    DB), this class can be used to download a bunch of FASTA files, or just one
    given its accession.

    Some **OLD** pre-defined lists are available from ENA. We refer to them as *virus*,
    *plasmid*, *phage*, *archaealvirus*, *archaea*, *bacteria*, *organelle*,
    *viroid*.

    .. warning:: the header of the FASTA files are changed to add the GI number
        instead of embl so th&at it can be used by our kraken builder class.
    """

    def __init__(self):
        """.. rubric:: constructor"""
        self.eutils = EUtilsTools()
        self.ena = ENA(cache=True)

    def ena_id_to_gi_number(self, identifiers):

        # Now, let us convert the ENA accession to NCBI GI number once for all.
        # We can fetch only at max 200 identifiers:
        logger.info("Fetching %s identifiers from NCBI" % len(identifiers))
        Nbaskets = int(math.ceil(len(identifiers) / 200.0))
        results = {}
        from easydev import split_into_chunks

        for chunk in split_into_chunks(identifiers, Nbaskets):
            result = self.eutils.accession_to_info(",".join(chunk))
            results.update(result)
        return results

    def download_fasta(self, filelist, output_dir=None):
        """Download a FASTA (or list of)

        :param filelist: a name to find on the ENA web server OR the
            name of an accession number or a file with accession numbers (1 column)

        .. warning:: The filename is named after the accession without .X number
            If there are several variant .1, .2 the later will be used. This
            should not happen if the list is properly defined.
        """
        try:
            # input may be a file
            with open(filelist, "r") as fin:
                self._identifiers = [x.strip() for x in fin.readlines()]
        except (FileNotFoundError, TypeError):
            try:
                # or a list of accession numbers
                self._identifiers = filelist + []
            except TypeError:
                # a just one accesion number
                self._identifiers = [filelist]

        self.results = self.ena_id_to_gi_number(self._identifiers)

        # do not use caching things this could be huge data sets.

        if output_dir is None:  # pragma: no cover
            output_dir = "."
        else:
            try:
                os.mkdir(output_dir)
            except FileExistsError:
                pass

        N = len(self._identifiers)
        pb = Progress(N)
        logger.info("Fetching all fasta from ENA")
        for i, identifier in enumerate(self._identifiers):
            # download data from ENA
            data = self.ena.get_data(identifier, "fasta")

            # Split header and Fasta
            try:
                header, others = data.split("\n", 1)
            except AttributeError:
                logger.warning(f"{identifier} not found on ENA website")
                continue

            # Source of failure:
            # - list and DB are not synchrone: e.g. some entries may be deleted
            if "suppressed" in header:  # pragma: no cover
                continue

            # Do not use try/except since when it fails, this is a real issue
            name = header.strip(">").split(" ")[0]
            db, id_, acc = name.split("|")

            header = self.add_gi_to_header(acc)

            # Save to local file
            # WARNINGS: extension is .fa because kraken-build expects .fa files
            filename = "%s_%s.fa" % (db, acc.split(".")[0])
            if output_dir:
                filename = output_dir + os.sep + filename

            with open(filename, "w") as fout:
                fout.write(header + "\n" + others)
            pb.animate(i + 1)

    def add_gi_to_header(self, acc):
        """Kraken will only accept the GI from NCBI so we need to convert
        the ENA accession to GI numbers"""

        # Accession may have a version .1, .2 hence this try/except first
        # without the version and then with the version.
        # Note also that some accession are different from an earlier version.
        # For instance, AF525933 is in the virus.txt list from ENA but
        # the new updated accession ois AH012103 showing that the list and DB
        # must not be fully synchronised.
        # http://www.ebi.ac.uk/ena/data/search?query=AF525933
        # In such case, the results attribute will be missing that accession,
        # which needs to be searched for specifically. We cannot now its name
        # before downloading the fasta.
        if acc in self.results.keys():
            res = self.results[acc]
            return f">{acc}|gi:{res['gi']}|taxid:{res['taxid']} {res['comment']}"
        else:  # pragma no cover
            logger.error(f"accession number not found on EUtils {acc}")
            raise KeyError


class NCBIDownload:
    def __init__(self):
        self.category = {
            "archaea",
            "bacteria",
            "fungi",
            "invertebrate",
            "mitochondrion",
            "other",
            "plant",
            "plasmid",
            "plastid",
            "protozoa",
            "vertebrate_mammalian",
            "vertebrate_other",
            "viral",
        }

        self.category_genomes = {
            "archaea",
            "bacteria",
            "fungi",
            "invertebrate",
            "plant",
            "protozoa",
            "vertebrate_mammalian",
            "vertebrate_other",
            "viral",
        }

    def download_ncbi_refseq_release(self, category, email="sequana@pasteur.fr"):
        """Download all files of type *fna* from ncbi FTP.

        ::

            kb = NCBIDownload()
            kb.download_ncbi_refseq_release("viral")

        """
        assert category in self.category, "Please use one of {}".format(self.category)

        ftp = ftplib.FTP("ftp.ncbi.nlm.nih.gov")
        ftp.login("anonymous", email)
        ftp.cwd(f"refseq/release/{category}")

        filenames = ftp.nlst()

        for filename in filenames:
            if "genomic.fna" in filename:
                logger.info(f"Downloading {filename}")
                ftp.retrbinary(f"RETR {filename}", open(filename, "wb").write)
        return [x for x in filenames if "genomic.fna" in x]

    def download_genomes_from_ncbi(
        self, category, email="sequana@pasteur.fr"
    ):  # pragma: no cover
        """This downloads all genomes on ncbi for a given category looking at
        their ftp. This could be highly redundant.

        """
        assert category in self.category

        ftp = ftplib.FTP("ftp.ncbi.nlm.nih.gov")
        ftp.login("anonymous", email)
        ftp.cwd(f"refseq/release/{category}")

        for filename in ftp.nlst():
            if "genomic.fna" in filename:
                logger.info(f"Downloading {filename}")
                ftp.retrbinary(f"RETR {filename}", open(filename, "wb").write)

    def download_assembly_report(self, category, output=None):
        assert category in self.category_genomes

        ftp = ftplib.FTP("ftp.ncbi.nlm.nih.gov")
        ftp.login("anonymous", "anonymous")
        ftp.cwd(f"genomes/refseq/{category}")

        filename = output if output else f"assembly_summary_{category}.txt"

        logger.info(f"Downloading {filename}")
        ftp.retrbinary(f"RETR assembly_summary.txt", open(filename, "wb").write)


class NCBITaxonReader:
    """This class will help in reading, handling, simplifying NCBI taxonomic DB

    When downloading NCBI taxonomy DB using e.g. Kraken, we end up with very
    large files. One is called names.dmp and the other nodes.dmp. They may be
    instrospected or simplified using this class

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
    and rank. Rank is species, genus, family, phylum, etc. Newest version of
    nodes.dmp hqs only 4 columns (taxid, parent taxid, rank ,a dash)

    ::

        from sequana.databases import NCBITaxonReader

        # The first time you may want to download the taxdump files
        n = NCBITaxonReader()
        n.download_taxdump()
        n.init("names.dmp", "nodes.dmp")

        # next time, you can read it directly
        n.NCBITaxonReader("names.dmp", "nodes.dmp")

    """

    ftp_url = "ftp.ncbi.nih.gov"

    def __init__(self, names=None, nodes=None):
        """.. rubric:: Constructor

        :param str names:  Defaults to "names.dmp".
        :param str nodes:  Defaults to "nodes.dmp".
        """
        if (names and os.path.exists(names)) and (nodes and os.path.exists(nodes)):
            self.init(names, nodes)
        else:
            logger.warning(
                "no input file provided or do not exist. Call download_taxdump() and init() manually"
            )

    def init(self, names, nodes):
        logger.info(f"Reading {names}")
        self.filename_names = names
        self.filename_nodes = nodes

        self.df_names = pd.read_csv(self.filename_names, sep="\t", header=None)
        # get rid of the pipes using ::2
        self.df_names = self.df_names.loc[:, ::2]
        self.df_names.rename(
            {0: "taxon", 2: "name", 4: "nothing", 6: "scname"}, axis=1, inplace=True
        )

        # This will provide a faster lookup table to search for scientific
        # names given a taxon. We can drop rows that are not scientific names
        # and set the taxons as index
        _subdf = self.df_names.query("'scientific name' in scname")
        self._subdf = _subdf.set_index("taxon")

        # Here, this is for general purpose (slower it we were to use
        # this for the get_scientic_name method
        self._group_name = self.df_names.groupby("taxon").groups

        logger.info(f"Reading {nodes}")

        self.df_nodes = pd.read_csv(self.filename_nodes, sep="\t", header=None)
        # get rid of the pipes using ::2
        self.df_nodes = self.df_nodes.loc[:, ::2]
        new_cols = ["taxon", "parent_taxon", "rank"] + list(
            range(0, len(self.df_nodes.columns) - 3)
        )
        self.df_nodes.columns = new_cols

        self._df_nodes_taxon = self.df_nodes.copy()
        self._df_nodes_taxon.set_index("taxon", inplace=True)

    def download_taxdump(self, outpath="."):  # pragma: no cover
        execute(
            f"wget {self.ftp_url}/pub/taxonomy/taxdump.tar.gz --directory-prefix {outpath}"
        )
        execute(f"tar xvfz {outpath}/taxdump.tar.gz -C {outpath}")

    def get_number_taxon(self):
        """Return number of unique taxon"""
        return len(self.df_names["taxon"].unique())

    def get_average_name_per_taxon(self):
        """Return number of rows/names per node/taxon"""
        from sequana.lazy import numpy as np

        return np.mean([len(values) for values in self._group_name.values()])

    def get_scientific_name(self, taxon):
        """Return scientific name of a given Taxon"""
        # Takes 2 minutes to scan all taxons
        return self._subdf.loc[taxon].values[0]

    def get_taxon_from_scientific_name(self, scname):
        """Return taxon corresponding to a scientific name

        return: unique taxon or first one found. If none found, returns None
        """
        res = self.df_names.query("@scname in name")["taxon"]
        return res

    def search(self, name):
        """Search names colum"""
        return self.df_names[self.df_names["name"].apply(lambda x: name in x)]

    def get_family(self, taxon):
        """Get all parent taxons"""
        taxons = [1]
        df = self._df_nodes_taxon

        final_links = (0, 1)
        while True:
            res = df.loc[taxon]
            taxons.append(taxon)
            taxon = res["parent_taxon"]
            # hopefully there is always a family link to 0 or 1
            if len(res) == 0 or taxon in final_links:
                break
        return taxons

    def filter_nodes_dmp_file(self, output="nodes_filtered.dmp", taxons=[]):
        """Save a subset of nodes.dmp given list of valid taxons

        :param str:  Defaults to "nodes_filtered.dmp".
        :param list taxons:
        """
        all_taxons = set()
        for i, taxon in enumerate(taxons):
            parents = self.get_family(taxon)
            all_taxons.update(parents)
        N = len(all_taxons)
        print(f"Saving {N} taxons")
        print(all_taxons)

        # In theory one line per taxon
        with open(self.filename_nodes, "r") as fin:
            with open(output, "w") as fout:
                for line in fin.readlines():
                    if int(line.split("\t", 1)[0]) in all_taxons:
                        fout.write(line)

    def filter_names_dmp_file(self, output="names_filtered.dmp", taxons=[]):
        """Save a subset of nodes.dmp given list of valid taxons

        :param str:  Defaults to "nodes_filtered.dmp".
        :param list taxons:
        """
        # here we want to save the family tree
        all_taxons = set()
        for i, taxon in enumerate(taxons):
            parents = self.get_family(taxon)
            all_taxons.update(parents)
        N = len(all_taxons)
        print(f"Saving {N} taxons")
        print(all_taxons)

        # we may have more entries than expecte because a same taxon is usesd for difference species
        with open(self.filename_names, "r") as fin:
            with open(output, "w") as fout:
                for line in fin.readlines():
                    if int(line.split("\t", 1)[0]) in all_taxons:
                        fout.write(line)
