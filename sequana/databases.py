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
"""
http://www.ebi.ac.uk/genomes/archaea.html

"""
import os
import glob
import math

try:
    # py2
    from urllib.request import urlopen
except:
    from urllib import urlopen


from bioservices import EUtils, ENA
from easydev import AttrDict


class EUtilsTools(object):

    def __init__(self):
        self.eutils = EUtils()

    def accession_to_info(self, ids):
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
    """

    Downloader facility to retrieve genome fasta files from ENA


    """
    def __init__(self):
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
            "macaca_fascicularis": ("macaca", "MacacaFascicularis")
        }

    def download_list(self):
        from easydev import execute
        for key, values in self._metadata.items():
            execute("wget http://www.ebi.ac.uk/genomes/%s" % values[0])

    def ena_id_to_gi_number(self, identifiers):
        # Now, let us convert the ENA accession to NCBI GI number once for all.
        # We can fetch only at max 200 identifiers:
        print("Fetching %s identifiers from NCBI" % len(identifiers))
        Nbaskets = math.ceil(len(identifiers)/200.)
        results = {}
        from easydev import split_into_chunks
        for chunk in split_into_chunks(identifiers, Nbaskets):
            result = self.eutils.accession_to_info(",".join(chunk))
            results.update(result)
        return results

    def download_fasta(self, filelist, output_dir=None, from_ena=True):
        """

        :param filelist: a name to find on the ENA web server OR the
            name of an accession number.


        .. warning:: The filename is named after the accession without .X number
            If there are several variant .1, .2 the later will be used. This
            should not happen if the list is properly defined. 
        """
        from bioservices import ENA
        if filelist.endswith(".txt") and from_ena is True:
            print("Downloading list from http://www.ebi.ac.uk/genomes/%s" % filelist)
            data = urlopen("http://www.ebi.ac.uk/genomes/%s" % filelist).readlines()
            identifiers = [x.strip().decode() for x in data]
        elif filelist == "macaca":
            identifiers = [ "CM001276", "CM001277", "CM001278", "CM001279",
                "CM001280", "CM001281", "CM001282", "CM001283", "CM001284", 
                "CM001285", "CM001286", "CM001287", "CM001288", "CM001289", 
                "CM001290", "CM001291","CM001292",  "CM001293", "CM001294", 
                "CM001295", "CM001296"]
        elif isinstance(filelist, list):
            pass
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

        from easydev import Progress
        N = len(identifiers)
        pb = Progress(N)
        print("Fetching all fasta from ENA")
        for i, identifier in enumerate(identifiers):
            filenames = glob.glob(output_dir + os.sep + "ENA_%s*" % identifier)

            if len(filenames) >= 1:
                pb.animate(i+1)
                # no need to fetch and save the data it looks like...
                continue

            # download data from ENA
            # http://www.ebi.ac.uk/ena/data/view/CM001276&display=fasta&download=fasta&filename=CM001276.fasta
            data = ena.get_data(identifier, "fasta")

            # Split header and Fasta
            header, others = data.decode().split("\n", 1)


            # Here, we will fetch info from NCBI every time. This is a pity but
            # no choice for now
            name  = header.strip(">").split(" ")[0]
            db, id_, acc = name.split("|")
            header = self.switch_header_to_gi(acc)
            # HEre we could have a try/except but better to fails since it
            # indicates in general that something is wrong e.g. unknown
            # identifier

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
        try:
            res = self.results[acc]
        except:
            res = self.results[acc.split(".")[0]]
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

        THis is the longest and takes about 20mins on a very good connection
        couple of hours otherwise.
        """
        self.download_fasta("bacteria.txt", "Bacteria")

    def download_accession(self, acc):
        self.download_fasta(acc, "Custom")

