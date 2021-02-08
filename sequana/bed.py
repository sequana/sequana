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
import colorlog
logger = colorlog.getLogger(__name__)



__all__ = ['BED']


class BED(object):
    """a structure to read and manipulate BED files (12-column file) 

    columns are defined as chromosome name, start and end, gene_name, score,
    strand, CDS start and end, blcok count, block sizes, block starts:

    """

    def __init__(self, filename):
        self.filename = filename

    def _get_line(self, line):
        try:
            if line.startswith('#') or  line.startswith('track') \
                or line.startswith('browser'):
                return {}
            else:
                fields = line.rstrip('\r\n').split()
                chrom_start = int(fields[1])
                return {"chrom": fields[0],
                        "start": chrom_start,
                        "end": int(fields[2]),
                        "gene_name": fields[3],
                        "score": fields[4],
                        "strand": fields[5],
                        "cds_start": int(fields[6]),
                        "cds_end": int(fields[7]),
                        "block_count": int(fields[9]),
                        "block_sizes": [int(x) for x in fields[10].strip(',').split(',')],
                        "block_starts": [int(x)+chrom_start for x in fields[11].strip(",").split(",")]}
        except Exception as err:
            print(err)
            print("Input bed must be 12-column] skipped line {}".format(line))
            return {}

    def __len__(self):
        count = 0
        with open(self.filename, "r") as fin:
            for line in fin:
                if line.startswith(('#', 'track', 'browser')): continue
                count += 1
        return count

    def get_exons(self):
        """Extract exon regions from input BED file.

        Uses the first (chromosome name), second (chromosome start), 11th and
        12th columns (exon start and size) of a 12-columns BED file.

        ::

            from sequana import sequana_data
            from sequana import BED
            b = BED(sequana_data("hg38_chr18.bed"))
            b.get_exons()

        """
        exons = []
        with open(self.filename, "r") as fin:
            for line in fin:
                if line.startswith('#'): continue
                if line.startswith('track'): continue
                if line.startswith('browser'): continue
                fields = line.rstrip('\r\n').split()
                assert len(fields) == 12 
                chrom_start = int(fields[1])
                chrom_name = fields[0]
                block_sizes =(int(x) for x in fields[10].strip(',').split(','))
                block_starts =  (int(x)+chrom_start for x in fields[11].strip(",").split(","))
                for start, size in zip(block_starts, block_sizes):
                    exons.append((chrom_name, start, start+size))
        return exons

    def get_transcript_ranges(self):
        """Extract transcript from input BED file."""
        with open(self.filename, "r") as fin:
            for line in fin:
                if line.startswith('#'): continue
                if line.startswith('track'): continue
                if line.startswith('browser'): continue
                fields = line.rstrip('\r\n').split()
                assert len(fields) == 12 
                chrom = fields[0]
                chrom_start = int(fields[1])
                chrom_end = int(fields[2])
                strand = fields[5]
                gene_name = fields[3]
                yield([chrom, chrom_start, chrom_end, strand,
                        "{}:{}:{}-{}".format(gene_name, chrom, chrom_start, chrom_end)])

    def get_CDS_exons(self):
        """Extract CDS from input BED file."""
        results = []
        with open(self.filename, "r") as fin:
            for line in fin:
                row = self._get_line(line)
                if row:
                    cds_start = row["cds_start"]
                    cds_end = row["cds_end"]
                    for start, size in zip( row['block_starts'], row['block_sizes'] ):
                        if (start + size) < cds_start: 
                            continue
                        if start > cds_end: 
                            continue
                        exon_start = max( start, cds_start )
                        exon_end = min( start+size, cds_end )
                        results.append([row['chrom'], exon_start, exon_end])
        return results



































