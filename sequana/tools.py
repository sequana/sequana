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
"""General tools"""
import os
import string
import glob
import json
import re
import gzip
import io
from collections import Counter

from sequana.lazy import pandas as pd
from sequana.lazy import numpy as np
from sequana import BAM

from pysam import FastxFile
from easydev import precision
from easydev.misc import cmd_exists
import subprocess

import colorlog
logger = colorlog.getLogger(__name__)



__all__ = ['StatsBAM2Mapped', 'bam_to_mapped_unmapped_fastq', "GZLineCounter"]


class DataContainer(dict):
    def __init__(self, wkdir="."):
        self.wkdir = wkdir

    def to_json(self, filename):
        json.dump(self.data, open(filename, "w"))

    def read_json(self, filename):
        return json.load(open(filename, "r"))

    def to_html(self):
        pass

try:
    _translate = string.maketrans('ACGTacgt', 'TGCAtgca')
    def bytes(x,dummy):
        return x
except:
    _translate = bytes.maketrans(b'ACGTacgt', b'TGCAtgca')


def reverse_complement(seq):
    # use hastag but for now, do it manually
    return seq.translate(_translate)[::-1]

def reverse(seq):
    return seq[::-1]


class StatsBAM2Mapped(DataContainer):

    def __init__(self, bamfile=None, wkdir=None, verbose=True):
        super(StatsBAM2Mapped, self).__init__(wkdir=wkdir)
        if bamfile.endswith(".bam"):
            self.data = bam_to_mapped_unmapped_fastq(bamfile, wkdir, verbose)
        elif bamfile.endswith(".json"):
            self.data = self.read_json(bamfile)

    def to_html(self):
        data = self.data

        # TODO hardcoded word phix here ?
        html = "Reads with Phix: %s %%<br>" % precision(data['contamination'], 3)

        # add HTML table
        if "R2_mapped" in data.keys():
            df = pd.DataFrame({
              'R1': [data['R1_mapped'], data['R1_unmapped']],
              'R2': [data['R2_mapped'], data['R2_unmapped']]})
        else:
            df = pd.DataFrame({
              'R1': [data['R1_mapped'], data['R1_unmapped']]})
        df.index = ['mapped', 'unmapped']

        html += "Unpaired: %s <br>" % data['unpaired']
        html += "duplicated: %s <br>" % data['duplicated']
        return html


def bam_to_mapped_unmapped_fastq(filename, output_directory=None, verbose=True):
    """Create mapped and unmapped fastq files from a BAM file

    :context: given a reference, one or two FastQ files are mapped onto the
        reference to generate a BAM file. This BAM file is a compressed version
        of a SAM file, which interpretation should be eased within this
        function.

    :param filename: input BAM file
    :param output_directory: where to save the mapped and unmapped files
    :return: dictionary with number of reads for each file (mapped/unmapped for
        R1/R2) as well as the mode (paired or not), the number of unpaired
        reads, and the number of duplicated reads. The unpaired reads should
        be zero (sanity check)

    Given a BAM file, create FASTQ with R1/R2 reads mapped and unmapped.
    In the paired-end case, 4 files are created.

    Note that this function is efficient in that it does not create intermediate
    files limiting IO in the process. As compared to standard tools such as 
    bedtools bamtofastq, it is 1.5 to 2X slower but it does create the mapped
    AND unmapped reads.

    :Details: Secondary alignment (flag 256) are dropped so as to remove any
        ambiguous alignments. The output dictionary stores "secondary" key to
        keep track of the total number of secondary reads that are dropped. If
        the flag is 256 and the read is unpaired, the key *unpaired* is also
        incremented.

        If the flag is not equal to 256, we first reverse complement reads that
        are tagged as *reverse* in the BAM file. Then, reads that are not paired or
        not "proper pair" (neither flag 4 nor flag 8) are ignored.

        If R1 is mapped **or** R2 is mapped then the reads are considered mapped. If
        both R1 and R2 are unmapped, then reads are unmapped.

    .. note:: about chimeric alignment: one is the representative and the other is
        the supplementary. This flag is not used in this function. Note also that
        chimeric alignment have same QNAME and flag 4 and 8

    .. note:: the contamination reported is basde on R1 only.

    .. todo:: comments are missing since there are not stored in the BAM file.


    .. note:: the mapped reads may not be synchronized because we include also
        the chimeric alignment (cf samtools documentation). However, 
        total reads = unmappeds reads + R1 mapped + R2 mapped - supplementary
        reads (those with flag 2048).
    """
    bam = BAM(filename)
    # figure out if this is paired or unpaired

    newname, ext = os.path.splitext(filename)

    import collections
    stats = collections.defaultdict(int)
    stats['R1_unmapped'] = 0
    stats['R1_mapped'] = 0

    # figure out where to save the file
    if output_directory is None:
        pass
    else:
        assert isinstance(filename, str)
        from sequana_pipetools.snaketools import FileFactory
        ff = FileFactory(filename)
        newname = output_directory + os.sep + ff.filenames[0]

    rt1 = "_R1_"
    rt2 = "_R2_"

    R1_mapped = open(newname + "{}.mapped.fastq".format(rt1), "wb")
    R1_unmapped = open(newname + "{}.unmapped.fastq".format(rt1), "wb")
    stats['duplicated'] = 0
    stats['unpaired'] = 0

    unpaired = 0

    # if paired, let open other files
    if bam.is_paired:
        stats['mode'] = "pe"
        stats['R2_unmapped'] = 0
        stats['R2_mapped'] = 0
        R2_mapped = open(newname + "{}.mapped.fastq".format(rt2), "wb")
        R2_unmapped = open(newname + "{}.unmapped.fastq".format(rt2), "wb")
    else:
        stats['mode'] = "se"

    # loop through the BAM (make sure it is rewinded)
    bam.reset()

    if verbose:
        from easydev import Progress
        pb = Progress(len(bam))

    for i, this in enumerate(bam):
        if this.flag & 256:
            # Unmapped reads are in the BAM file but have no valid assigned
            # position (N.B., they may have an assigned position, but it should be ignored).
            # It's typically the case that a number of reads can't be aligned, due to things
            # like sequencing errors, imperfect matches between the DNA sequenced and the
            # reference, random e. coli or other contamination, etc..
            # A secondary alignment occurs when a given read could align reasonably well to
            # more than one place. One of the possible reported alignments is termed "primary"
            # and the others will be marked as "secondary".
            stats['secondary'] += 1
            if this.is_paired is False:
                stats['unpaired'] += 1
        else:

            # quick hack
            if this.is_read1:
                suffix = b"/1"
            else:
                suffix = b"/2"

            # in pysam, seq is a string and qual a bytes....
            if this.is_reverse is True:
                txt = b"@" + bytes(this.qname, "utf-8") + suffix + b"\n"
                revcomp = reverse_complement(this.seq)
                txt += bytes(revcomp, "utf-8") + b"\n"
                txt += b"+\n"
                txt += bytes(this.qual[::-1], 'utf-8') + b"\n"
            else:
                txt = b"@" + bytes(this.qname, "utf-8") + suffix + b"\n"
                txt += bytes(this.seq, "utf-8") + b"\n"
                txt += b"+\n"
                txt += bytes(this.qual,"utf-8") + b"\n"

            # Here, we must be careful as to keep the pairs. So if R1 is mapped
            # but R2 is unmapped (or the inverse), then the pair is mapped
            if this.is_read1:
                if this.is_unmapped and this.mate_is_unmapped:
                    R1_unmapped.write(txt)
                    stats['R1_unmapped'] += 1
                else:
                    R1_mapped.write(txt)
                    stats['R1_mapped'] += 1
            elif this.is_read2:
                if this.is_unmapped and this.mate_is_unmapped:
                    R2_unmapped.write(txt)
                    stats['R2_unmapped'] += 1
                else:
                    R2_mapped.write(txt)
                    stats['R2_mapped'] += 1
            else:
                # This should be a single read
                #assert self.is_paired is False
                stats['unpaired'] += 1
                if this.is_unmapped:
                    R1_unmapped.write(txt)
                    stats['R1_unmapped'] += 1
                else:
                    R1_mapped.write(txt)
                    stats['R1_mapped'] += 1

            if this.is_duplicate:
                stats['duplicated'] += 1

        if verbose:
            pb.animate(i+1)

    if bam.is_paired:
        R2_mapped.close()
        R2_unmapped.close()

    logger.info("\nNumber of entries in the BAM: %s" % str(i+1))

    R1_mapped.close()
    R1_unmapped.close()

    _x = stats['R1_mapped']
    _y = stats['R1_unmapped']
    stats["contamination"] = _x / float(_x + _y) * 100

    return stats


def bam_get_paired_distance(filename):
    """Return distance between 2 mated-reads

    :return: list of tuples where each tuple contains the position start,
        position end of the paired-end reads that were mapped + the mode.
        mode =1 means fragment is reversed. mode = 2 means mate is reversed.
        mode = 3 means none are reversed.

    ::

        distances = bam_get_paired_distance(bamfile)
        hist([x[1]-x[0] for x in distances])

    .. warning:: experimental
    """

    b = BAM(filename)
    distances = []

    for fragment in b:
        if fragment.is_unmapped is False and fragment.mate_is_unmapped is False \
            and fragment.is_read1:

            # get the mate:
            mate = next(b)


            if fragment.is_reverse:
                position2 = fragment.reference_end
                position1 = mate.reference_start
                mode = 1
            elif mate.is_reverse:
                position1 = fragment.reference_start
                position2 = mate.reference_end
                mode = 2
            else: # if both are not reversed, what does that mean.
                # On Hm2, this is the case for 4 pairs out of 1622
                # This seems to be a special case for fragment ends exactly
                # at the end of the reference and mate starts exactly at
                # the beginnin with a length less than 100
                print(fragment.reference_start, fragment.reference_end)
                print(mate.reference_start, mate.reference_end)
                position1 = -1
                position2 = -1
                mode = 3

            distances.append((position1, position2, mode))

    return distances


def _base_content(filename, window_size, letters, circular=False):
    # DOC: see gc_content
    fasta = FastxFile(filename)
    checker = set(letters)
    chrom_gc_content = dict()
    for chrom in fasta:
        mid = int(window_size / 2)
        # Create gc_content array
        gc_content = np.empty(len(chrom.sequence))
        gc_content[:] = np.nan
        if circular:
            chrom.sequence = (chrom.sequence[-mid:] + chrom.sequence +
                    chrom.sequence[:mid])
            # Does not shift index of array
            mid = 0
        # Count first window content
        counter = Counter(chrom.sequence[0:window_size])
        gc_count = 0
        for letter in letters:
            gc_count += counter[letter]

        gc_content[mid] = gc_count
        for i in range(1, len(chrom.sequence) - window_size + 1):
            if chrom.sequence[i - 1] in checker:
                gc_count -= 1
            if chrom.sequence[i + window_size - 1] in checker:
                gc_count += 1
            gc_content[i + mid] = gc_count
        chrom_gc_content[chrom.name] = gc_content / window_size
    return chrom_gc_content


def gc_content(filename, window_size, circular=False, 
        letters=['G', 'C', 'c', 'g']):
    """Return GC content for the different sequences found in a FASTA file

    :param filename: fasta formated file
    :param window_size: window length used to compute GC content
    :param circular: set to True if sequences are circular.
    :return: dictionary with keys as fasta names and values as GC content vecor

    .. todo:: case when the genome is not circular -> Add NaN at start and stop of
        the np.arange()

    """
    return _base_content(filename, window_size, letters, circular=circular)


def genbank_features_parser(input_filename):
    """ Return dictionary with features contains inside a genbank file.

    :param str input_filename: genbank formated file
    """
    logger.warning("deprecated. please use GenBank.genbank_features_parser instead")
    new_feature = {}
    records = {}
    feature_list = []
    feature_field = False

    with open(input_filename, "r") as fp:
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




class GZLineCounter(object):
    """Fast GZipped line counter

    Uses zcat if possible, otherwise gzip library (twice as slow).

    .. doctest::

        >>> from sequana import sequana_data
        >>> from sequana.misc import GZLineCounter
        >>> gz = GZLineCounter(sequana_data("test.fastq.gz"))
        >>> len(gz)
        100

    """
    def __init__(self, filename):
        self.filename = filename
        if cmd_exists("zcat"):
            self.use_zcat = True
        else:
            self.use_zcat = False

    def __len__(self):
        if self.use_zcat:
            return self._use_zcat()
        else:
            return self._use_gzip()

    def _use_zcat(self):
        i = 0
        p = subprocess.Popen(["zcat", self.filename], stdout=subprocess.PIPE)
        for line in p.stdout:
            i +=1
        return i

    """twice as fast as the other method below but large memory print
    def _use_gzip(self):
        with gzip.open(self.filename) as gz_file:
            data = gz_file.read()
        return len(data.splitlines())
    """
    def _use_gzip(self):
        i = 0
        with gzip.open(self.filename) as gz_file:
            with io.BufferedReader(gz_file) as f:
                for line in f:
                    i += 1
        return i


def entropy(x):
    letters = set(x)
    from pylab import log
    pi = [x.count(l)/float(len(x)) for l in letters]
    pi = [x for x in pi if x!=0]
    return -sum(pi*log(pi))




class PairedFastQ(object):

    def __init__(self, fq1, fq2):
        self.fq1 = fq1
        self.fq2 = fq2

    def is_synchronised(self):
        from sequana import FastQ
        N = 0
        for a, b in zip(FastQ(self.fq1), FastQ(self.fq2)):
            a = a['identifier'].decode()
            b = b['identifier'].decode()
            a = a.split()[0]
            b = b.split()[0]

            if a.endswith("/1"):
                id1 = a.rsplit("/1")[0]
            elif a.endswith("/2"):
                id1 = a.rsplit("/2")[0]
            else:
                id1 = a
            if b.endswith("/1"):
                id2 = b.rsplit("/1")[0]
            elif b.endswith("/2"):
                id2 = b.rsplit("/2")[0]
            else:
                id2 = b

            if id1 != id2:
                print("%s differs from %s" % (id1, id2))
                print(a)
                print(b)
                return False
            N += 1
        print(N)
        return True
