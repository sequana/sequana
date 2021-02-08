# -*- coding: utf-8 -*-
#
#  This file is part of Sequana software
#
#  Copyright (c) 2016 - Sequana Development Team
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
"""Utilities for GTF files"""
from collections import defaultdict
import sys

from sequana.fastq import FastQ

import colorlog
logger = colorlog.getLogger(__name__)




class GTFFixer():
    """Some GTF have syntax issues. This class should fix some of them

    With exons, we want to have unique IDs.

    We append the transcrpt version and ID to the exon ID
    """
    def __init__(self, filename):

        self.filename = filename

    def fix(self, output):
        features = defaultdict(int)
        logger.info("Scanning file")
        count = 0
        exon_ids = []
        with open(self.filename, "r") as fin:
            with open(output, "w") as fout:
                line = fin.readline()
                count += 1
                while line:
                    if line.startswith("#") or len(line.strip()) == 0:
                        fout.write(line)
                    else:
                        entries = line.split(maxsplit=8)
                        features[entries[2]] += 1
                        try:
                            # sanity check that this is correctly written for all
                            # features
                            annotations = dict([x.strip().split(maxsplit=1) for x in entries[8].split(";") if len(x.strip())])
                        except:
                            logger.warning("warning line {}. could not parse annotations correctly".format(count))

                        # if exon feature, we want to make sure the ID are unique by
                        # appending the transcript_version if not already done
                        if entries[2] == "exon" and "." not in annotations["exon_id"]:
                            logger.info("Fixing exon ID line {} for exon ID {}".format(count, annotations["exon_id"]))
                            exon_id = "{}.{}.{}".format(annotations['exon_id'].rstrip('"'), annotations['transcript_version'].strip('"'), annotations['transcript_id'].lstrip('"'))
                            newline = "\t".join(entries[0:8])+"\t"
                            # we do not resuse the dictionaru 'annotations' just to
                            # keep the order
                            for entry in entries[8].split(";"):
                                if len(entry.strip()):
                                    x, y = entry.split(maxsplit=1)
                                    if x.strip() != "exon_id":
                                        newline += "{} {};".format(x,y)
                                    else:
                                        newline += "exon_id {};".format(exon_id)
                            newline +=";\n"
                            line = newline

                        fout.write(line)
                    line = fin.readline()
                    count += 1
                    if count % 100000 == 0: print("scanned {} entries".format(count))
            return features


