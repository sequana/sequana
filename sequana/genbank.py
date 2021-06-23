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
import re
from sequana.fasta import FastA

import colorlog
logger = colorlog.getLogger(__name__)



__all__ = ["GenBank"]


# TODO: we should factorise gff and genbank in a parent class (Annotation)
class GenBank():
    """
    ::

        gg = GenBank()
        gg.get_types()

    """
    def __init__(self, filename):
        self.filename = filename

    def get_types(self):

        records = self.genbank_features_parser()
        _types = set()
        for contig in records.keys():
            for feature in records[contig]:
                _type = feature['type']
                _types.add(_type)
        return sorted(_types)

    def extract_fasta(self, fastafile, features=['rRNA']):
        types = self.get_types()
        for feature in features:
            if feature not in types:
                raise ValueError("{} not found".format(feature))

        # fasta may have several contig/chromosome names
        # the gene bank should be compatible !!
        fasta = FastA(fastafile)
        contig_names = fasta.get_lengths_as_dict()

        # most of the times, the version is not in the gbk
        contig_names = [x.split(".")[0] for x in contig_names]

        # then we read the features from the genbank
        records = self.genbank_features_parser()
        contig_names_gbk = list(records.keys())

        # FIXME FastA is not very efficient for eukaryotes but is enough for now

        output = ""
        for name in records.keys():
            if name not in contig_names:
                logger.warning("{} contig from genbank not found in fasta".format(name))
                continue
            index = contig_names.index(name)
            sequence = fasta.sequences[index]

            for item in records[name]:
                if item['type'] in features:
                    start, end = item['gene_start'], item['gene_end']
                    try:
                        info = item['product']
                        output += ">{}_{}_{}_{} {}\n".format(name, item['type'], 
                                                             start,end, info)
                    except:
                        output += ">{}_{}_{}_{} {}\n".format(name, item['type'], start, end)
                    output+= "{}\n".format(sequence[start:end])
        return output

    def genbank_features_parser(self):
        """ Return dictionary with features contains inside a genbank file.

        :param str input_filename: genbank formated file
        """
        new_feature = {}
        records = {}
        feature_list = []
        feature_field = False

        with open(self.filename, "r") as fp:
            for line in fp:
                # pass header and sequence fields
                if not feature_field:
                    # get contig/chrom name
                    if line.startswith("LOCUS"):
                        name = line.split()[1]
                    elif line.startswith("FEATURE"):
                        feature_field = True
                else:
                    # if feature field is finished
                    if line.startswith("ORIGIN"):
                        feature_field = False
                        records[name] = feature_list
                        feature_list = []
                        new_feature = []
                        continue

                    # if there are a word in qualifier indent (feature type)
                    # maybe we need to infer the size of indentation ???
                    if line[0:20].split():
                        if new_feature:
                            feature_list.append(new_feature)
                        split_line = line.split()
                        t = split_line[0]
                        # Handle :
                        #complement(order(1449596..1449640,1449647..1450684,
                        #1450695..1450700))
                        positions = split_line[1]
                        if positions[0].isalpha():
                            while not line[:-1].endswith(")"):
                                line = next(fp)
                                positions += line
                        pos = [int(n) for n in re.findall(r"\d+", positions)]
                        # Handle complement(join(3773333..3774355,3774357..3774431))
                        start = pos[0]
                        end = pos[-1]
                        strand = "-" if split_line[1].startswith("c") else "+"
                        new_feature = {"type": t, "gene_start": start,
                                "gene_end": end, "strand": strand}
    
                    # recover qualifier bound with feature
                    else:
                        quali_line = line.strip().replace('"', '')
                        if quali_line.startswith("/") and "=" in quali_line:
                            qualifier = quali_line.split("=")
                            key = qualifier[0][1:]
                            new_feature[key] = qualifier[1]
                        else:
                            if key == "translation":
                                new_feature[key] += quali_line
                            else:
                                new_feature[key] += " " + quali_line
        return records


