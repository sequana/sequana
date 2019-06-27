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
"""Utilities to access to online FASTA, taxon, lineage ..."""
import os
import glob
import math

try:
    # py2
    from urllib.request import urlopen
except:
    from urllib import urlopen

from easydev import AttrDict, execute, Progress

from sequana import logger
logger.name = __name__


class EUtilsTools(object):
    """Simple wrapper around EUtils to fetch basic informatino about an accession number


    ::

        >>> from sequana.databases import EUtilsTools
        >>> et.accession_to_info("K01711.1")
        {'K01711.1': {'accession': '331784',
          'comment': 'Measles virus (strain Edmonston), complete genome',
          'gi': '331784',
          'identifier': 'gi|331784|gb|K01711.1|MEANPCG[331784]',
          'taxid': '11234'}}


    """
    def __init__(self):
        from bioservices import EUtils
        self.eutils = EUtils()

    def accession_to_info(self, ids):
        """An accession or list of them returns list of dictionaries"""
        res = self.eutils.EFetch(db="nuccore",id=ids,
                rettype="docsum", retmode="dict")

        res = res['eSummaryResult']['DocSum']

        # if one id provided, it will be a dict, otherwise a list of dicts
        try:
            res[0]
        except:
            res = [res]

        # now we can loop over all identifiers
        records = {}
        accessions = [x.strip() for x in ids.split(',')]

        for i, entry in enumerate(res):
            # first, save the acc number
            accession = entry['Id']
            # then various info
            items = entry['Item']
            identifier = [x for x in items if x['@Name'] == "Extra"][0]['#text']
            if "||" in identifier:
                # strip content after ||
                identifier = identifier.split("||")[0]

            title = [x for x in items if x['@Name'] == "Title"][0]['#text']
            taxid = [x for x in items if x['@Name'] == "TaxId"][0]['#text']
            gi = [x for x in items if x['@Name'] == "Gi"][0]['#text']
            record = {
                "taxid": taxid,
                'accession': accession,
                "identifier": identifier,
                'gi': gi,
                'comment': title
                }

            records[accessions[i]] = AttrDict(**record)
        return records


class ENADownload(object):
    """Downloader to retrieve genome fasta files from ENA amongst other things

    In order to facilitate the download of FASTA files (e.g. to build a Kraken
    DB), this class can be used to download a bunch of FASTA files, or just one
    given its accession. 

    Pre-defined lists are available from ENA. We refer to them as *virus*,
    *plasmid*, *phage*, *archaealvirus*, *archaea*, *bacteria*, *organelle*,
    *viroid*. In addition we have predefined lists within Sequana. For now,
    there is one named macaca_fascicularis.

    .. warning:: the header of the FASTA files are changed to add the GI number
        instead of embl.
    """
    def __init__(self):
        """.. rubric:: constructor"""
        self.convert_enaacc_to_gi = True
        self.eutils = EUtilsTools()

        # In the tuple, the first element is either a ebi txt files to obtain
        # via wget or locally, and the second item is the output directory
        self._metadata = {
            'virus': ("virus.txt", "Viruses"),
            'plasmid': ("plasmid.txt", "Plasmids"),
            'phage': ("phage.txt", "Phages"),
            'archaealvirus': ("archaealvirus.txt", "ArchaeaViruses"),
            'archaea': ("archaea.txt", "Archaea"),
            'bacteria': ("bacteria.txt", "Bacteria"),
            'organelle': ("organelle.txt", "Organelle"),
            'viroid': ("viroid.txt", "Viroid"),
            "macaca_fascicularis": ("macaca", "MacacaFascicularis"),
            "mus_musculus": ("mus_musculus", "MusMusculus"),
            "worms": ("worms", "Worms")
        }

    def download_list(self):
        """Download all standard lists of accession numbers from ENA"""
        for key, values in self._metadata.items():
            execute("wget -q -t 3 http://www.ebi.ac.uk/genomes/%s -O %s"  
                % (values[0], values[0]))

    def ena_id_to_gi_number(self, identifiers):

        # Now, let us convert the ENA accession to NCBI GI number once for all.
        # We can fetch only at max 200 identifiers:
        logger.info("Fetching %s identifiers from NCBI" % len(identifiers))
        Nbaskets = int(math.ceil(len(identifiers)/200.))
        results = {}
        from easydev import split_into_chunks
        for chunk in split_into_chunks(identifiers, Nbaskets):
            result = self.eutils.accession_to_info(",".join(chunk))
            results.update(result)
        return results

    def download_fasta(self, filelist, output_dir=None, from_ena=True):
        """Download a FASTA (or list of)

        :param filelist: a name to find on the ENA web server OR the
            name of an accession number.

        .. warning:: The filename is named after the accession without .X number
            If there are several variant .1, .2 the later will be used. This
            should not happen if the list is properly defined. 
        """
        from bioservices import ENA
        if filelist.endswith(".txt") and os.path.exists(filelist) is False:
            logger.info("Downloading list from http://www.ebi.ac.uk/genomes/%s" % filelist)
            data = urlopen("http://www.ebi.ac.uk/genomes/%s" % filelist).readlines()
            identifiers = [x.strip().decode() for x in data]
        elif filelist == "macaca":
            identifiers = [ "CM001276", "CM001277", "CM001278", "CM001279",
                "CM001280", "CM001281", "CM001282", "CM001283", "CM001284",
                "CM001285", "CM001286", "CM001287", "CM001288", "CM001289",
                "CM001290", "CM001291","CM001292",  "CM001293", "CM001294",
                "CM001295", "CM001296"]
        elif filelist == "mus_musculus": #19 +x+y chromosomes + 5 mitochondrion
            # could also add strain C57BL.
            identifiers = ["AY172335", "CM000209", "CM000210", "CM000211"
                "CM000212", "CM000213", "CM000214", "CM000215", "CM000216"
                "CM000217", "CM000218", "CM000219", "CM000220", "CM000221"
                "CM000222", "CM000223", "CM000224", "CM000225", "CM000226"
                "CM000227", "CM000228", "CM000229", "CM000225", "CM000226"
                "EF108342", "AB042432", "AY675564", "DQ874614"]
        elif filelist == "worms": # Caernorhabditis briggsae and elegans
            identifiers = ["AC186293", "FR847112", "FR847113", "FR847114",
                "FR847118", "FR847121", "FR847123", "BX284601", "BX284602",
                "BX284603", "BX284604", "BX284605", "BX284606"]
        elif isinstance(filelist, str) and filelist in self._metadata.keys():
            name = self._metadata[filelist][0]
            logger.info("Downloading list from http://www.ebi.ac.uk/genomes/%s" % name)
            data = urlopen("http://www.ebi.ac.uk/genomes/%s" % name).readlines()
            identifiers = [x.strip().decode() for x in data]
        elif isinstance(filelist, list):
            identifiers = filelist[:]
        elif isinstance(filelist, str):
            # could be a single identifier or a filename (assuming a single
            # column)
            if os.path.exists(filelist):
                identifiers = [x for x in open(filelist).read().split()]
                identifiers = [x.strip() for x in identifiers]
            else:
                identifiers = [filelist]
        self._identifiers = identifiers

        self.results = self.ena_id_to_gi_number(identifiers)

        # do not use caching things this could be huge data sets.
        ena = ENA()

        if output_dir is None:
            output_dir = "."
        else:
            try: os.mkdir(output_dir)
            except:pass

        N = len(identifiers)
        pb = Progress(N)
        logger.info("Fetching all fasta from ENA")
        for i, identifier in enumerate(identifiers):
            filenames = glob.glob(output_dir + os.sep + "ENA_%s*" % identifier)

            if len(filenames) >= 1:
                pb.animate(i+1)
                # no need to fetch and save the data it looks like...
                continue

            # download data from ENA
            data = ena.get_data(identifier, "fasta")

            # Split header and Fasta
            header, others = data.decode().split("\n", 1)

            # Source of failure:
            # - list and DB are not synchrone: e.g. some entries may be deleted
            if "suppressed" in header:
                continue
            if ">" not in header:
                continue

            # Do not use try/except since when it fails, this is a real issue
            name  = header.strip(">").split(" ")[0]
            db, id_, acc = name.split("|")

            try:
                header = self.switch_header_to_gi(acc)
            except:
                logger.error("Failed for this entry:") 
                logger.error(identifier)
                logger.error(header)
                logger.error(name)
                continue

            # Save to local file
            # WARNINGS: extension is .fa because kraken-build expects .fa files
            filename = "%s_%s.fa" % (db, acc.split(".")[0])
            if output_dir:
                filename = output_dir + os.sep + filename

            with open(filename, "w") as fout:
                fout.write(header+"\n"+others)
            pb.animate(i+1)

    def switch_header_to_gi(self, acc):
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
        else:
            try:
                res = self.results[acc.split(".")[0]]
            except:
                logger.warning("\nUnknown accession (%s). May be an updated version. Checking..." % acc)
                res = self.ena_id_to_gi_number([acc])
                self.results.update(res)
                res = res[acc]
                logger.info('Found %s using GI number' % acc)
        return ">"+res['identifier']+" " + res['comment']

    def download_viroid(self):
        self.download_fasta(*self._metadata['viroid'])

    def download_organelle(self):
        self.download_fasta(*self._metadata['organelle'])

    def download_viruses(self):
        self.download_fasta(*self._metadata['virus'])

    def download_plasmids(self):
        self.download_fasta("plasmid.txt", "Plasmids")

    def download_phage(self):
        self.download_fasta("phage.txt", "Phages")

    def download_archaealvirus(self):
        self.download_fasta("archaealvirus.txt", "ArchaeaViruses")

    def download_archaea(self):
        self.download_fasta(*self._metadata["archaea"])

    def download_macaca(self):
        self.download_fasta(*self._metadata["macaca_fascicularis"])

    def download_bacteria(self):
        """ organisms (may 2016)

        .. note:: this download method is the longest to end. It took about 20mins on 
            a good connection.
        """
        self.download_fasta("bacteria.txt", "Bacteria")

    def download_accession(self, acc, output="Custom"):
        """Download a specific FASTA file given its ENA accession number """
        self.download_fasta(acc, output)

