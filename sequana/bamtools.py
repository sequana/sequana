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
"""Tools to manipulate BAM/SAM files

.. autosummary::

    Alignment
    BAM
    CRAM
    MultiBAM
    SAM
    SAMFlags
    SAMBAMbase 

.. note:: BAM being the compressed version of SAM files, we do not
    implement any functionalities related to SAM files. We strongly encourage
    developers to convert their SAM to BAM.

"""
import os
import json
import math

from collections import Counter, OrderedDict, defaultdict

from bx.intervals.intersection import Interval, IntervalTree
from bx.bitset import BinnedBitSet 
from bx.bitset_builders import binned_bitsets_from_list
from sequana.lazy import pandas as pd
from sequana.lazy import numpy as np
from sequana.lazy import pylab

from sequana.bed import BED
from sequana.cigar import fetch_intron, fetch_exon

import pysam
from sequana import jsontool

import colorlog
logger = colorlog.getLogger(__name__)


"""
#http://www.acgt.me/blog/2014/12/16/understanding-mapq-scores-in-sam-files-does-37-42#
#http://biofinysics.blogspot.fr/2014/05/how-does-bowtie2-assign-mapq-scores.html
#https://gitlab.univ-nantes.fr/a-slide/ContaVect/blob/9a411abfa720064c205c5f6c811afdfea206ed12/pyDNA/pySamTools/Bam.py

Interesting commands::

    samtools flagstat contaminant.bam
"""

# SAMBAMbase is for the doc
__all__ = ['BAM','Alignment', 'SAMFlags', "CS", "SAM", "CRAM", "SAMBAMbase"]


# simple decorator to rewind the BAM file
# here we use an underscore to not overlap with the reset() method
# of the AlignmentFile class
from functools import wraps
def _reset(f):
    @wraps(f)
    def wrapper(*args, **kargs):
        args[0].reset()
        return f(*args, **kargs)
    return wrapper


# There are lots of trouble with inheriting from pysam.AlignmentFile
# First, you cannot use super(). indeed, it works for py27 but not 
# with py35 probably a missing __init__  or __new__ in
# AlignmentFile class. See e.g., stackoverflow/questions/
# 26653401/typeerror-object-takes-no-parameters-but-only-in-python-3
# So we should call the __init__.

# Second problem. The parent method reset() works for BAM but not for SAM 
# files. Indeed, the reset() method rewind the file but with SAM, it seems to 
# be confused with the header somehow. So we need to overload the reset()
# to call the __init__ again but this is not appropriate to call the constructor
# again. It may have side effects. 

# So, we decided to store the SAM/BAM as an attribute. 

# Other known issues with pysam.AlignmentFile:
# - in PY3, there is an attribute **format** that tells us if the input is BAM
# or SAM but this does not work in PY2 where the attribute is missing...
# - we cannot use the method is_bam trustfully.
# - Using AlignmentFile one must provide the proper mode 'e.g. rb will set
# is_bam to True while 'r' only will set it to False. However, it seems
# there is not sanity check inside the file.


def is_bam(filename, *args):
    """Return True if input file looks like a BAM file"""
    f = pysam.AlignmentFile(filename, mode="r", *args)
    return f.is_bam


def is_sam(filename, *args):
    """Return True if input file looks like a SAM file"""
    f = pysam.AlignmentFile(filename, mode="r", *args)
    return f.is_sam


def is_cram(filename, *args):
    """Return True if input file looks like a CRAM file"""
    f = pysam.AlignmentFile(filename, mode="r", *args)
    return f.is_cram


class SAMBAMbase():
    """Base class for SAM/BAM/CRAM data sets


    We provide a few test files in Sequana, which can be retrieved with
    sequana_data:

    .. doctest::

        >>> from sequana import BAM, sequana_data
        >>> b = BAM(sequana_data("test.bam"))
        >>> len(b)
        1000
        >>> from sequana import CRAM
        >>> b = CRAM(sequana_data("test_measles.cram"))
        >>> len(b)
        60

    """
    # The mode rb means read-only (r) and that (b) for binary the format
    # So BAM or SAM can be read in theory.
    def __init__(self, filename, mode="r", *args):
        self._filename = filename
        self._mode = mode
        self._args = args
        self._summary = None
        self._sorted = None

        # Save the length so that second time we need it, it is already
        # computed.
        self._N = None
        self.reset()

        # we can accumulate the length of each contig/chromosome as an alias
        bam = pysam.AlignmentFile(self._filename, mode=self._mode, *self._args)
        hd = bam.header
        self.lengths = {r: l for r,l in zip(hd.references, hd.lengths)}

    def reset(self):
        try:
            self._data.close()
        except:
            pass
        self._data = pysam.AlignmentFile(self._filename,
            mode=self._mode, *self._args)

    @_reset
    def get_read_names(self):
        """Return the reads' names"""
        names = [this.qname for this in self._data]
        return names

    @_reset
    def __len__(self):
        if self._N is None:
            logger.warning("Scanning the BAM. Please wait")
            self._N = sum(1 for _ in self._data)
            self.reset()
        return self._N

    @_reset
    def _get_insert_size_data(self, max_entries=100000):
        count = 0
        data = []
        for a in self:
            count += 1
            if a.is_paired and a.tlen!=0:
                data.append(a.tlen)
            if count>max_entries:
                break
        return data

    @_reset
    def get_estimate_insert_size(self, max_entries=100000, upper_bound=1000,
        lower_bound=-1000):
        """

        .. image:: _static/insert_size.png

        Here we show that about 3000 alignments are enough to get a good
        estimate of the insert size.

        .. plot::

            from sequana import *
            from pylab import linspace, plot, grid, xlabel, ylabel

            b = BAM(sequana_data("measles.fa.sorted.bam"))
            X = linspace(100,3000,1000)
            plot(X, [b.get_estimate_insert_size(x) for x in X])
            grid()
            xlabel("Number of alignements used")
            ylabel("Estimated insert size")

        """
        # Get an estimate of the read length
        data = self._get_insert_size_data(max_entries=max_entries)
        if len(data) == 0:
            return 0
        data = [abs(x) for x in data if x>=lower_bound and x<=upper_bound]
        return pylab.mean([abs(x) for x in data])

    @_reset
    def get_df(self, max_align=-1):
        flags = []
        starts = []
        ends = []
        mapqs = []
        refnames = []
        querynames = []
        querylengths = []
        cigar = []
        for i, a in enumerate(self._data):
            flags.append(a.flag)
            starts.append(a.reference_start)
            ends.append(a.reference_end)
            mapqs.append(a.mapq)
            cigar.append(a.cigarstring)
            try:
                refnames.append(a.reference_name)
            except:
                refnames.append(-1)
            querynames.append(a.query_name)
            querylengths.append(a.query_length)
            if max_align!=-1 and i>max_align:
                break
        df = pd.DataFrame({
                "flag": flags,
                "rstart": starts,
                "rend": ends,
                "mapqs": mapqs,
                "rname": refnames,
                "qname": querynames,
                "qlen": querylengths
            })
        return df

    @_reset
    def get_df_concordance(self, max_align=-1):
        """This methods returns a dataframe with Insert, Deletion, Match,
        Substitution, read length, concordance (see below for a definition)


        Be aware that the SAM or BAM file must be created using minimap2 and the
        --cs option to store the CIGAR in a new CS format, which also contains
        the information about substitution. Other mapper are also handled (e.g.
        bwa) but the substitution are solely based on the NM tag if it exists.

        alignment that have no CS tag or CIGAR are ignored.


        """
        from sequana import Cigar
        count = 0
        I, D, M, L, mapq, flags, NM, rnames = [], [], [], [], [], [], [], []
        S = []
        for i, a in enumerate(self._data):
            # tags and cigar populated  if there is a match
            # if we use --cs cigar is not populated so we can only look at tags
            # tags can be an empty list
            if a.tags is None or len(a.tags) == 0:
                continue
            count += 1
            mapq.append(a.mapq)
            L.append(a.qlen)
            try:
                NM.append([x[1] for x in a.tags if x[0] == "NM"][0])
            except:
                #FIXME why -1 and not 0
                NM.append(-1)

            flags.append(a.flag)
            rnames.append(a.reference_name)

            if 'cs' in dict(a.tags):
                cs = CS(dict(a.tags)['cs'])
                S.append(cs['S'])
                I.append(cs['I'])
                D.append(cs['D'])
                M.append(cs['M'])
            elif a.cigarstring:
                cigar = Cigar(a.cigarstring).as_dict()
                I.append(cigar["I"])
                D.append(cigar['D'])
                M.append(cigar['M'])
                S.append(None)  # no info about substitutions in the cigar
            else:
                I.append(0)
                D.append(0)
                M.append(0)
                S.append(0)

            if max_align>0 and count == max_align:
                break

            if count % 10000 == 0:
                logger.debug("Read {} alignments".format(count))

        I = np.array(I)
        D = np.array(D)
        M = np.array(M)
        NM = np.array(NM)
        rnames = np.array(rnames)

        try:
            S = np.array(S)
            denom = S + I + D + M
            C = 1 - (I + D + S) / denom
            logger.info("computed Concordance based on minimap2 --cs option")
        except:
            logger.info("computed Concordance based on standard CIGAR information using INDEL and NM tag")
            computed_S = NM - D - I
            C = 1 - (I + D + computed_S)/(computed_S + I + D + M)

        df = pd.DataFrame([C, L, I, D, M, mapq, flags, NM, S, rnames])
        df = df.T
        df.columns = ["concordance", 'length', "I", "D", "M", "mapq", "flags", "NM", "mismatch", "rname"]
        return df

    def __iter__(self):
        return self

    def __next__(self):
        return next(self._data)

    @_reset
    def infer_strandness(self, reference_bed, max_entries, mapq=30):
        """
        :param reference_bed: a BED file (12-columns with 
            columns 1,2,3,6 used) or GFF file (column 1, 3, 
            4, 5, 6 are used
        :param mapq: ignore alignment with mapq below 30.
        :param max_entries: can be long. max_entries restrict the estimate

        Strandness of transcript is determined from annotation while
        strandness of reads is determined from alignments.

        For non strand-specific RNA-seq data, strandness of reads and strandness
        of transcript are independent.

        For strand-specific RNA-seq data, strandness of reads is determined by
        strandness of transcripts.
        

        This functions returns a list of 4 values. First one indicates whether
        data is paired or not. Second and third one are ratio of reads explained
        by two types of strandness of reads vs transcripts. Last values are
        fractions of reads that could not be explained. The values 2 and 3 tell
        you whether this is a strand-specificit dataset. 

        If similar, it is no strand-specific.  If the first value is close to 1
        while the other is close to 0, this is a strand-specific dataset

        """
        count = 0
        p_strandness = defaultdict(int)
        s_strandness = defaultdict(int)

        gene_ranges = {}
        if reference_bed.endswith(".bed"):
            with open(reference_bed, "r") as fin:
                for line in fin:
                    if line.startswith(('#','track','browser')):
                        continue
                    # Parse fields from gene tabls
                    fields = line.split()
                    chrom_name = fields[0]
                    start = int( fields[1] )
                    end = int( fields[2] )
                    #gene_name = fields[3]
                    strand = fields[5]

                    if chrom_name not in gene_ranges:
                        gene_ranges[chrom_name] = IntervalTree()
                    gene_ranges[chrom_name].insert(start,end, strand)
        elif reference_bed.endswith(".gff"):
            with open(reference_bed, "r") as fin:
                count = 0
                for line in fin:
                    count +=1
                    if line.startswith(('#')):
                        continue
                    fields = line.split()
                    if len(fields)<6:
                        logger.warning("invalid format on line {}: {}".format(count, line))
                    if fields[2] != "gene":
                        continue
                    chrom_name = fields[0]
                    start = int( fields[3] )
                    end = int( fields[4] )
                    #gene_name = fields[3]
                    strand = fields[6]
                    if chrom_name not in gene_ranges:
                        gene_ranges[chrom_name] = IntervalTree()
                    gene_ranges[chrom_name].insert(start,end, strand)

        self.gene_ranges = gene_ranges

        count = 0
        for aln in self:
            if count >= max_entries:
                break
            if aln.is_qcfail:
                continue
            if aln.is_duplicate:
                continue
            if aln.is_secondary:
                continue
            if aln.is_unmapped:
                continue
            if aln.mapq < mapq:
                continue

            chrom_name = self._data.get_reference_name(aln.tid)

            if aln.is_paired:
                if aln.is_read1:
                    read_id = '1'
                if aln.is_read2:
                    read_id = '2'
                if aln.is_reverse:
                    map_strand = '-'
                else:
                    map_strand = '+'
                readStart = aln.pos
                readEnd = readStart + aln.qlen
                if chrom_name in gene_ranges:
                    tmp = set(gene_ranges[chrom_name].find(readStart, readEnd))
                    if len(tmp) == 0: 
                        continue
                    strand_from_gene = ':'.join(tmp)
                    p_strandness[read_id + map_strand + strand_from_gene] += 1
                    count += 1
            else:
                if aln.is_reverse:
                    map_strand = '-'
                else:
                    map_strand = '+'
                readStart = aln.pos
                readEnd = readStart + aln.qlen
                if chrom_name in gene_ranges:
                    tmp = set(gene_ranges[chrom_name].find(readStart,readEnd))
                    if len(tmp) == 0: continue
                    strand_from_gene = ':'.join(tmp)
                    s_strandness[map_strand + strand_from_gene] += 1
                    count += 1

        self.s_strandness = s_strandness
        self.p_strandness = p_strandness
        # Fraction of reads failed to determine: 0.0189
        # Fraction of reads explained by "1++,1--,2+-,2-+": 0.6315
        # Fraction of reads explained by "1+-,1-+,2++,2--": 0.3497
        # 
        if len(p_strandness) >0 and len(s_strandness) ==0 :
            protocol = "Paired-end"
            spec1= (p_strandness['1++'] + p_strandness['1--'] 
                    + p_strandness['2+-'] + p_strandness['2-+'])/float(sum(p_strandness.values()))
            spec2= (p_strandness['1+-'] + p_strandness['1-+'] 
                    + p_strandness['2++'] + p_strandness['2--'])/float(sum(p_strandness.values()))
            other = 1 - spec1 - spec2

        elif len(s_strandness) > 0 and len(p_strandness) == 0 :
            protocol = "Singled-end"
            spec1 = (s_strandness['++'] + s_strandness['--']) / float(sum(s_strandness.values()))
            spec2 = (s_strandness['+-'] + s_strandness['-+']) / float(sum(s_strandness.values()))
            other = 1-spec1-spec2
        else:
            protocol = "Mixture"
            spec1 = "NA"
            spec2 = "NA"
            other = "NA"
        return [protocol, spec1, spec2, other]

    @_reset
    def mRNA_inner_distance(self, refbed, low_bound=-250, up_bound=250,
            step=5, sample_size=1000000, q_cut=30):

        """Estimate the inner distance of mRNA pair end fragment.

        ::

            from sequana import BAM, sequana_data
            b = BAM(sequana_data("test_hg38_chr18.bam"))
            df = b.mRNA_inner_distance(sequana_data("hg38_chr18.bed"))


        """
        #This code was inspired from the RSeQC code v2.6.4 and adapted for
        #sequana simplifying the code and using pandas to store results. This is
        #limited by memory but more convenient.

        inner_distance_bitsets = BinnedBitSet()
        tmp = BinnedBitSet()
        tmp.set_range(0,0)
        pair_num=0

        inner_distances = []
        read_names = []
        descriptions = []

        logger.info("Get exon regions from ")

        bed_obj = BED(refbed)
 
        ref_exons = [[exn[0].upper(), exn[1], exn[2]] for exn in bed_obj.get_exons()]

        exon_bitsets = binned_bitsets_from_list(ref_exons)

        transcript_ranges = {}
        for i_chr, i_st, i_end, i_strand, i_name in bed_obj.get_transcript_ranges():
            i_chr = i_chr.upper()
            if i_chr not in transcript_ranges:
                transcript_ranges[i_chr] = IntervalTree()
            else:
                transcript_ranges[i_chr].add_interval(Interval(i_st, i_end, value=i_name))

        try:
            while(1):
                if pair_num >= sample_size:
                    break
                splice_intron_size=0
                aligned_read = next(self)

                if aligned_read.is_qcfail:continue           #skip low quality
                if aligned_read.is_duplicate:continue        #skip duplicate read
                if aligned_read.is_secondary:continue        #skip non primary hit
                if aligned_read.is_unmapped:continue         #skip unmap read
                if not aligned_read.is_paired: continue      #skip single map read
                if aligned_read.mate_is_unmapped:continue    #
                if aligned_read.mapq < q_cut:continue

                read1_len = aligned_read.qlen
                read1_start = aligned_read.pos
                read2_start = aligned_read.mpos        #0-based, not included
                if read2_start < read1_start:
                    continue                           #because BAM file is sorted, mate_read is already processed if its coordinate is smaller
                if  read2_start == read1_start and aligned_read.is_read1:
                    inner_distance = 0
                    #inner_distances.append(0)
                    continue

                pair_num +=1

                # check if reads were mapped to diff chromsomes
                R_read1_ref = self._data.get_reference_name(aligned_read.tid)
                R_read2_ref = self._data.get_reference_name(aligned_read.rnext)
                if R_read1_ref != R_read2_ref:
                    inner_distances.append("NA")
                    read_names.append(aligned_read.qname)
                    descriptions.append("sameChrom=No")
                    continue

                chrom = self._data.get_reference_name(aligned_read.tid).upper()
                intron_blocks = fetch_intron(chrom, read1_start, aligned_read.cigar)
                for intron in intron_blocks:
                    splice_intron_size += intron[2] - intron[1]
                read1_end = read1_start + read1_len + splice_intron_size

                if read2_start >= read1_end:
                    inner_distance = read2_start - read1_end
                else:
                    exon_positions = []
                    exon_blocks = fetch_exon(chrom, read1_start,aligned_read.cigar)
                    for ex in exon_blocks:
                        for i in range(ex[1]+1,ex[2]+1):
                            exon_positions.append(i)
                    inner_distance = -len([i for i in exon_positions if i > read2_start and i <= read1_end])

                read1_gene_names = set()    #read1_end
                try:
                    for gene in transcript_ranges[chrom].find(read1_end-1, read1_end):
                        #gene: Interval(0, 10, value=a)
                        read1_gene_names.add(gene.value)                        
                except:
                    pass

                read2_gene_names = set()    #read2_start
                try:
                    for gene in transcript_ranges[chrom].find(read2_start, read2_start +1):
                        #gene: Interval(0, 10, value=a)
                        read2_gene_names.add(gene.value)
                except:
                    pass

                if len(read1_gene_names.intersection(read2_gene_names)) == 0:
                    # no common gene
                    #reads mapped to different gene
                    inner_distances.append(inner_distance)
                    read_names.append(aligned_read.qname)
                    descriptions.append("sameTranscript=No,dist=genomic")
                    continue

                if inner_distance > 0:
                    if chrom in exon_bitsets:
                        size =0
                        inner_distance_bitsets.set_range(read1_end, read2_start-read1_end)
                        inner_distance_bitsets.iand(exon_bitsets[chrom])
                        end=0
                        while 1:
                            start = inner_distance_bitsets.next_set( end )
                            if start == inner_distance_bitsets.size: break
                            end = inner_distance_bitsets.next_clear( start )
                            size += (end - start)
                        #clear BinnedBitSet
                        inner_distance_bitsets.iand(tmp)   

                        if size == inner_distance:
                            inner_distances.append(size)
                            read_names.append(aligned_read.qname)
                            descriptions.append("sameTranscript=Yes,sameExon=Yes,dist=mRNA")
                        elif size > 0 and size < inner_distance:
                            inner_distances.append(size)
                            read_names.append(aligned_read.qname)
                            descriptions.append("sameTranscript=Yes,sameExon=No,dist=mRNA")
                        elif size <= 0:
                            inner_distances.append(inner_distance)
                            read_names.append(aligned_read.qname)
                            descriptions.append("sameTranscript=Yes,nonExonic=Yes,dist=genomic")
                    else:
                        inner_distances.append(inner_distance)
                        read_names.append(aligned_read.qname)
                        descriptions.append("unknownChromosome,dist=genomic")
                else:
                    inner_distances.append(inner_distance)
                    read_names.append(aligned_read.qname)
                    descriptions.append("readPairOverlap")

        except StopIteration:
            pass

        print("Total read pairs  used " + str(pair_num))
        if pair_num==0:
            raise ValueError("Cannot find paired reads")

        df = pd.DataFrame(
            {"read_names": read_names,
             "val": inner_distances,
             "desc": descriptions})
        df = df[["read_names", "val", "desc"]]
        df.query("val>=@low_bound and val<=@up_bound").val.hist(bins=50)
        mu = df.query("val>=@low_bound and val<=@up_bound").val.mean()
        print("mean insert size: {}".format(mu))
        pylab.title("Mean inner distance={}".format(round(pylab.mean(mu), 2) ))
        pylab.axvline(mu, color="r", ls="--", lw=2)
        return df, mu

    # properties
    @_reset
    def _get_paired(self):
        return next(self).is_paired
    is_paired = property(_get_paired)

    @_reset
    def _get_is_sorted(self):
        if self._sorted:
            return self._sorted
        pos = next(self._data).pos
        for this in self._data:
            if this.is_unmapped is True:
                continue
            if this.pos < pos:
                self._sorted = False
                return False
            pos = this.pos
        self._sorted = True
        return self._sorted
    is_sorted = property(_get_is_sorted, doc="return True if the BAM is sorted")

    def get_samtools_stats_as_df(self):
        """Return a dictionary with full stats about the BAM/SAM file

        The index of the dataframe contains the flags. The column contains
        the counts.

        ::

            >>> from sequana import BAM, sequana_data
            >>> b = BAM(sequana_data("test.bam"))
            >>> df = b.get_samtools_stats_as_df()
            >>> df.query("description=='average quality'")
            36.9

        .. note:: uses samtools behind the scene
        """
        from easydev import shellcmd
        res = shellcmd("samtools stats %s" % self._filename)
        res = res.decode('utf-8')

        # First, we can extract all data that statrts with SN
        # The format is
        #
        # SN name: value #comment
        #
        # separators are \t tabulation
        #
        # so we split with the : character, remove the starting SN\t characters
        # remove comments and ignore other \t characters. We should end up with
        # only 2 columns; names/values

        # extra all relevnt lines starting with SN
        data = [x for x in res.split("\n") if x.startswith('SN')]

        # remove comments
        data = [x.split('#')[0][3:] for x in data]
        names = [x.split(":")[0] for x in data]
        values = [x.split(":")[1].strip() for x in data]
        df = pd.DataFrame({"description": names, "count": values })
        df = df[['description', 'count']]
        df.sort_values(by='count', inplace=True)
        return df

    def _count_item(self, d, item, n=1):
        if item in d.keys():
            d[item] += n
        else:
            d[item] = n

    def _get_summary(self):
        """Count flags/mapq/read length in one pass."""
        if self._summary is not None:
            return self._summary

        mapq_dict = {}
        read_length_dict = {}
        flag_dict = {}
        mean_qualities = []
        count = 0
        for read in self:
            self._count_item(mapq_dict, read.mapq)
            self._count_item(flag_dict, read.flag)
            if read.is_unmapped is False:
                self._count_item(read_length_dict, read.reference_length)
            try:
                mean_qualities.append(pylab.mean(read.query_qualities))
            except TypeError:
                mean_qualities.append(-1)
            count += 1
            if count % 100000 ==0:
                print(count)
        # FIXME do we need the try/except if so, add Exception
        try:
            mq = pylab.mean(mean_qualities)
        except:
            mq = 0
        self._summary = {"mapq": mapq_dict,
                         "read_length": read_length_dict,
                         "flags": flag_dict,
                         "mean_quality": pylab.mean(mean_qualities)
                         }
        return self._summary
    summary = property(_get_summary)

    def _get_read_length(self):
        X = sorted(self.summary['read_length'].keys())
        Y = [self.summary['read_length'][k] for k in X]
        return X, Y

    def plot_insert_size(self, max_entries=100000, bins=100, upper_bound=1000,
        lower_bound=-1000):
        """

        This gives an idea of the insert size without taking into account any
        intronic gap. The mode should give a good idea of the insert size
        though.


        .. plot::

            from sequana import *
            from pylab import linspace, plot, grid, xlabel, ylabel

            b = BAM(sequana_data("measles.fa.sorted.bam"))
            b.plot_insert_size()

        """
        data = self._get_insert_size_data(max_entries=max_entries)
        if len(data) == 0:
            return 0
        data = [x for x in data if x>=lower_bound and x<=upper_bound]
        M = pylab.mean([abs(x) for x in data])
        pylab.hist(data, bins=bins)
        return M

    def plot_read_length(self):
        """Plot occurences of aligned read lengths

        .. plot::
            :include-source:

            from sequana import sequana_data, BAM
            b = BAM(sequana_data("test.bam"))
            b.plot_read_length()

        """
        X, Y = self._get_read_length()
        pylab.plot(X, Y,
            label="min length:{}; max length:{}".format(min(X), max(X)))
        pylab.grid()
        pylab.xlabel("Read length", fontsize=16)
        pylab.legend()

    def get_stats(self):
        """Return basic stats about the reads

        :return: dictionary with basic stats:

            - total_reads : number reads ignoring supplementaty and secondary
              reads
            - mapped_reads : number of mapped reads
            - unmapped_reads : number of unmapped
            - mapped_proper_pair : R1 and R2 mapped face to face
            - reads_duplicated: number of reads duplicated

        .. warning:: works only for BAM files. Use :meth:`get_samtools_stats_as_df`
            for SAM files.

        """

        """#See samtools stats
        # 1526795 + 0 in total (QC-passed reads + QC-failed reads)
        13 + 0 secondary
        0 + 0 supplementary
        0 + 0 duplicates
        3010 + 0 mapped (0.20% : N/A)
        1526782 + 0 paired in sequencing
        763391 + 0 read1
        763391 + 0 read2
        2700 + 0 properly paired (0.18% : N/A)
        2976 + 0 with itself and mate mapped
        21 + 0 singletons (0.00% : N/A)
        0 + 0 with mate mapped to a different chr
        0 + 0 with mate mapped to a different chr (mapQ>=5)
        """
        d = {}

        samflags_count = self.get_samflags_count()

        # all reads - (supplementary alignmnt + secondary alignmnt)
        d['total_reads'] = len(self) - (samflags_count[256] +
                                        samflags_count[2048])
        # all reads - (unmapped + supplementary alignmnt + secondary alignmnt)
        d['mapped_reads'] = d['total_reads'] - samflags_count[4]
        d['unmapped_reads'] = samflags_count[4]
        d['mapped_proper_pair'] = samflags_count[2]
        d['reads_duplicated'] = samflags_count[1024]
        d['secondary_reads'] = samflags_count[256]
        return d

    @_reset
    def get_stats_full(self, mapq=30, max_entries=-1):
        # On a bam, this takes about 7minutes
        # while calling samtools directly takes a little bit more than 1 minute.
        # but then we can re-use this independently of samtools (not pysam though)

        average_quality = 0
        average_length = 0
        bases_mapped = 0
        bases_mapped_cigar = 0  # more precise according to samtools
        forward = 0
        reads_duplicated = 0
        insert_size_sum = 0
        insert_size_sum_square = 0
        non_splice = 0
        mapq0 = 0
        mismatches = 0
        multiple_hit = 0
        qc_fail = 0
        reverse = 0
        reads_paired = 0
        unique_hit = 0
        unmapped = 0
        read1 = 0
        read2 = 0
        splice = 0
        secondary = 0
        pair_diff_chrom = 0
        proper_pair = 0
        total_aln = 0
        total_length = 0
        total_r1_length = 0
        total_r2_length = 0

        for aln in self:
            total_aln += 1
            average_quality += sum(aln.query_qualities) / len(aln.query_qualities)

            if aln.is_paired:
                reads_paired += 1

            if aln.is_unmapped:
                unmapped +=1

            if aln.is_qcfail: 
                qc_fail +=1
                continue
            if aln.is_duplicate:
                reads_duplicated += 1
                continue
            if aln.is_secondary:
                secondary += 1
                continue

            total_length += aln.rlen

            # fixme not really a multiple hit in 100% of cases
            if aln.mapq < mapq:
                multiple_hit += 1
                continue
            else:
                unique_hit += 1
                if aln.is_read1:
                    read1 += 1
                    total_r1_length += aln.rlen
                if aln.is_read2:
                    read2 += 1
                    total_r2_length += aln.rlen
                if aln.is_reverse:
                    reverse += 1
                else:
                    forward +=1

            if aln.mapq == 0:
                mapq += 0

            # average length of all MAPPED reads to divide by read1+read2 that
            # mapped
            average_length += aln.rlen
            A = abs(aln.tlen)

            bases_mapped_cigar += sum([a[1] for a in aln.cigar if a[0] in [0,1]])
            bases_mapped += aln.rlen

            introns = fetch_intron("dummy", aln.pos, aln.cigar)
            if len(introns) == 0:
                non_splice += 1
            elif len(introns)>=1:
                splice += 1

            if aln.is_proper_pair:
                proper_pair += 1
                read1_ref = self._data.get_reference_name(aln.tid)
                read2_ref = self._data.get_reference_name(aln.rnext)
                if read1_ref != read2_ref:
                    pair_diff_chrom += 1
                if aln.pos != 0:
                    # to avoid effect of circular genome
                    insert_size_sum += A
                    insert_size_sum_square += A * A

            # If NM is provided, not always the case though
            try:
                mismatches += aln.get_tag('NM')
            except:
                pass

            if total_aln % 100000 == 0:
                print(total_aln)
                if max_entries != -1 and total_aln >= max_entries:
                    break
            if max_entries != -1 and total_aln >= max_entries:
                break

        results = {
            "average_quality": average_quality / total_aln,
            "average_length": average_length / (read1+read2),
            "bases mapped (cigar)":  bases_mapped_cigar,
            "bases mapped ":  bases_mapped,
            "forward": forward,
            "is sorted": self.is_sorted,
            "unmapped": unmapped,
            "mismatches": mismatches,
            "multiple_hit": multiple_hit,
            "non_splice": non_splice,
            "pair_diff_chrom": pair_diff_chrom,
            "proper_pair": proper_pair,
            "qc_fail": qc_fail,
            "read1": read1,
            "read2": read2,
            "reads_duplicated": reads_duplicated,
            "reads_mapq0": mapq0,
            "reads_mapped": read1 + read2,
            "reads_paired": reads_paired,
            # In theory, the next one is paired-end technology bit set + both mates mapped
            "reads_mapped_and_paired": reads_paired - secondary,
            "raw_total_sequences": total_aln - secondary,
            "reverse": reverse,
            "non_primary_alignements": secondary,
            "secondary": secondary,
            "splice": splice,
            "total_alignments": total_aln,
            "total_first_fragment_length": total_r1_length,
            "total_last_fragment_length": total_r2_length,
            "total_length": total_length,  # ignoring secondary, qcfail to agree with samtools
            "unique_hit": unique_hit,
            "error_rate": mismatches / bases_mapped_cigar
            }
        assert results['forward'] + results['reverse']

        if proper_pair> 0:
            results["insert_size_average"] =  insert_size_sum  / proper_pair
            results["insert_size_std"] =  math.sqrt(
                insert_size_sum_square/(proper_pair) - (insert_size_sum/(proper_pair))**2)
        if reads_paired>0:
            results["percentage_properly_paired"] = 100*proper_pair / reads_paired
        else:
            results["percentage_properly_paired"] = 0
        #else:
        #    results["insert_size_std"] = None

        return results

    """
SN	insert size average:	4775.2
SN	insert size standard deviation:	3817.6
SN	inward oriented pairs:	169
SN	outward oriented pairs:	229
SN	pairs with other orientation:	3
SN	pairs on different chromosomes:	0


 'bases mapped ': 70050,
 'bases mapped (cigar)': 65641,
 'forward': 458,
 'insert_size_average': 3040185.341991342,
 'multiple_hit': 2,
 'total_length': 72347,


"""

# Note that on large files, differences can be large. Here on a human genome, we
# got: 
# 'insert_size_average': 2406.5781705775985,
# 'insert_size_std': 11633.214727733652,
#
#and with samtools:
#            insert size average        1219.9
# insert size standard deviation        2249.1
# In both case, this is wrong. This was RNA data sets and using mRNA_insert_size
# code, we get about 72 bases as expected. One should ignore values that are too
# large e.g. below and above -250/250

    @_reset
    def get_flags_as_df(self):
        """Returns decomposed flags as a dataframe

        .. doctest::

            >>> from sequana import BAM, sequana_data
            >>> b = BAM(sequana_data('test.bam'))
            >>> df = b.get_flags_as_df()
            >>> df.sum()
            0          0
            1       1000
            2        484
            4          2
            8          2
            16       499
            32       500
            64       477
            128      523
            256       64
            512        0
            1024       0
            2048       0
            dtype: int64

        .. seealso:: :class:`SAMFlags` for meaning of each flag
        """
        flags = [s.flag for s in self]
        data = [(this, [flag&this for flag in flags])
            for this in (0, 1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048)]
        df = pd.DataFrame(dict(data))

        # special case of flag 0 has to be handled separetely. Indeed 0 & 0 is 0
        # If flag is zero, we store 1, otherwise 0
        df[0] = [1 if x==0 else 0 for x in flags]

        df = df > 0
        return df

    def plot_bar_flags(self, logy=True, fontsize=16, filename=None):
        """Plot an histogram of the flags contained in the BAM

        .. plot::
            :include-source:

            from sequana import BAM, sequana_data
            b = BAM(sequana_data('test.bam', "testing"))
            b.plot_bar_flags()

        .. seealso:: :class:`SAMFlags` for meaning of each flag
        """
        df = self.get_flags_as_df()
        df = df.sum()
        pylab.clf()
        if logy is True:
            barplot = df.plot(kind='bar', logy=logy, grid=True)
        else:
            barplot = df.plot(kind='bar', grid=True)
        pylab.xlabel("flags", fontsize=fontsize)
        pylab.ylabel("count", fontsize=fontsize)
        pylab.tight_layout()
        if filename:
            pylab.savefig(filename)
        return barplot

    @_reset
    def to_fastq(self, filename):
        """Export the BAM to FastQ format

        .. todo:: comments from original reads are not in the BAM so will be missing

        Method 1 (bedtools)::

            bedtools bamtofastq -i JB409847.bam  -fq test1.fastq

        Method2 (samtools)::

            samtools bam2fq JB409847.bam > test2.fastq

        Method3 (sequana)::

            from sequana import BAM
            BAM(filename)
            BAM.to_fastq("test3.fastq")

        Note that the samtools method removes duplicated reads so the output is
        not identical to method 1 or 3.

        """
        with open(filename, "w") as fh:
            for i, this in enumerate(self):
                # FIXME what about comments. not stored in the BAM 
                read = this.qname
                read += this.seq + "\n"
                read += "+\n"
                read += this.qual + "\n"
                #if i != self.N-1:
                #    read += "\n"
                fh.write(read)

    @_reset
    def get_mapq_as_df(self, max_entries=-1):
        """Return dataframe with mapq for each read"""
        if max_entries != -1:
            data = [next(self).mapq for x in range(max_entries)]
        else:
            data = [this.mapq for this in self]
        df = pd.DataFrame({'mapq': data})
        return df

    @_reset
    def get_mapped_read_length(self):
        """Return dataframe with read length for each read


        .. plot::

            from pylab import hist
            from sequana import sequana_data, BAM
            b = BAM(sequana_data("test.bam"))
            hist(b.get_mapped_read_length())

        """
        read_length = [read.reference_length for read in self
                       if read.is_unmapped is False]
        return read_length

    def get_samflags_count(self):
        """ Count how many reads have each flag of SAM format.


        :return: dictionary with keys as SAM flags
        """
        samflags = (1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048)
        samflags_count = dict.fromkeys(samflags, 0)
        for flag, count in self.summary["flags"].items():
            for samflag in samflags:
                if flag & samflag != 0:
                    samflags_count[samflag] += count
        return samflags_count

    def plot_bar_mapq(self, fontsize=16, filename=None, ):
        """Plots bar plots of the MAPQ (quality) of alignments

            .. plot::
                :include-source:

                from sequana import BAM, sequana_data
                b = BAM(sequana_data('test.bam', "testing"))
                b.plot_bar_mapq()

        """
        df = self.get_mapq_as_df()
        df.plot(kind='hist', bins=range(0,df.max().values[0]+1), legend=False,
            grid=True, logy=True)
        pylab.xlabel("MAPQ", fontsize=fontsize)
        pylab.ylabel("Count", fontsize=fontsize)
        pylab.tight_layout()
        if filename:
            pylab.savefig(filename)

    def bam_analysis_to_json(self, filename):
        """ Create a json file with information related to the bam file.

        This includes some metrics (see :meth:`get_stats`; eg MAPQ),
        combination of flags, SAM flags, counters about the read length.
        """
        d = {}
        d["module"] = "bam_analysis"
        d["metrics"] = self.get_stats()
        d["combo_flag"] = self.summary["flags"]
        d["samflags"] = self.get_samflags_count()
        d["read_length"] = self.summary["read_length"]
        with open(filename, "w") as fp:
            json.dump(d, fp, indent=True, sort_keys=True)

    @_reset
    def get_gc_content(self):
        """Return GC content for all reads (mapped or not)

        .. seealso:: :meth:`plot_gc_content`

        """
        data = [(f.seq.count("C") + f.seq.count('G')) / len(f.seq)*100. for f in self if f.seq]
        return data

    @_reset
    def get_length_count(self):
        """Return counter of all fragment lengths"""
        import collections
        data = [this.rlen for this in self]
        return collections.Counter(data)

    def plot_gc_content(self, fontsize=16, ec="k", bins=100):
        """plot GC content histogram

        :params bins: a value for the number of bins or an array (with a copy()
            method)
        :param ec: add black contour on the bars

        .. plot::
            :include-source:

            from sequana import BAM, sequana_data
            b = BAM(sequana_data('test.bam'))
            b.plot_gc_content()

        """
        data = self.get_gc_content()
        try:
            X = np.linspace(0, 100, bins)
        except:
            X = bins.copy()

        pylab.hist(data, X, density=True, ec=ec)
        pylab.grid(True)
        mu = pylab.mean(data)
        sigma = pylab.std(data)

        X = pylab.linspace(X.min(), X.max(), 100)

        from sequana.misc import normpdf

        pylab.plot(X, normpdf(X, mu, sigma), lw=2, color="r", ls="--")
        pylab.xlabel("GC content", fontsize=16)

    def _get_qualities(self, max_sample=500000):
        qualities = []
        for i, record in enumerate(self):
            if i < max_sample:
                #quality = [ord(x) -33 for x in record.qual]
                quality = record.query_qualities
                qualities.append(quality)
            else:
                break
        return qualities

    @_reset
    def boxplot_qualities(self, max_sample=500000):
        """Same as in :class:`sequana.fastq.FastQC`

        """
        qualities = self._get_qualities(max_sample)
        df = pd.DataFrame([x for x in qualities if x])
        from sequana.viz.boxplot import Boxplot
        bx = Boxplot(df)
        try:
            bx.plot(ax=ax)
        except:
            bx.plot()

    #FIXME: why not a property ? Same comments for coverage attribute
    def _set_alignments(self):
        # this scans the alignments once for all
        self.alignments = [this for this in self]

    @_reset
    def _set_coverage(self):
        try:
            self.alignments
        except AttributeError:
            self._set_alignments()
        ref_start = defaultdict(list)
        ref_end = defaultdict(list)
        
        for aln in self.alignments:
            # Of course, we must have a valid reference name and start/end position not set to None
            if aln.flag not in [4] and aln.rname != -1 and aln.reference_end and aln.reference_start:
                ref_start[aln.rname].append(aln.reference_start)
                ref_end[aln.rname].append(aln.reference_end)

        self.coverage = {}

        for rname in ref_start.keys():
            print(rname)
            print(None in ref_end[rname])
            N = max(ref_end[rname])
            self.coverage[rname] = np.zeros(N)
            for x, y in zip(ref_start[rname], ref_end[rname]):
                if y and x>=0 and y>=0:
                    self.coverage[rname][x:y] += 1
                else:
                    pass

    @_reset
    def _set_indels(self):
        try:
            self.alignments
        except:
            self._set_alignments()

        self.insertions = []
        self.deletions = []
        for this in self.alignments:
            if this.cigarstring:
                if "I" in this.cigarstring:
                    self.insertions.extend([x[1] for x in this.cigartuples if x[0] == 1])
                if "D" in this.cigarstring:
                    self.deletions.extend([x[1] for x in this.cigartuples if x[0] == 2])

    def plot_coverage(self, chrom=None):
        """Please use :class:`GenomeCov` for more sophisticated
        tools to plot the genome coverage

        .. plot::
            :include-source:

            from sequana import sequana_data, BAM
            b = BAM(sequana_data("measles.fa.sorted.bam"))
            b.plot_coverage()

        """
        try:
            self.coverage
        except AttributeError:
            self._set_coverage()

        if chrom is None and len(self.coverage.keys()) == 1:
            chrom = list(self.coverage.keys())[0]

        pylab.plot(self.coverage[chrom])
        pylab.xlabel("Coverage")

    def hist_coverage(self, chrom=None, bins=100):
        """

        .. plot::
            :include-source:

            from sequana import sequana_data, BAM
            b = BAM(sequana_data("measles.fa.sorted.bam"))
            b.hist_coverage()
        """
        try:
            self.coverage
        except AttributeError:
            self._set_coverage()

        if chrom is None and len(self.coverage.keys()) == 1:
            chrom = list(self.coverage.keys())[0]

        pylab.hist(self.coverage[chrom], bins=bins)
        pylab.xlabel("Coverage")
        pylab.ylabel("Number of mapped bases")
        pylab.grid()

    @_reset
    def plot_indel_dist(self, fontsize=16):
        """Plot indel count (+ ratio)

        :Return: list of insertions, deletions and ratio insertion/deletion for
            different length starting at 1

        .. plot::
            :include-source:

            from sequana import sequana_data, BAM
            b = BAM(sequana_data("measles.fa.sorted.bam"))
            b.plot_indel_dist()

        What you see on this figure is the presence of 10 insertions of length
        1, 1 insertion of length 2 and 3 deletions of length 1


        # Note that in samtools, several insertions or deletions in a single
        alignment are ignored and only the first one seems to be reported. For
        instance 10M1I10M1I stored only 1 insertion in its report; Same comment
        for deletions.

        .. todo:: speed up and handle long reads cases more effitiently by 
            storing INDELS as histograms rather than lists
        """
        try:
            self.insertions
        except:
            self._set_indels()

        if len(self.insertions) ==0 or len(self.deletions) == 0:
            raise ValueError("No deletions or insertions found")

        N = max(max(Counter(self.deletions)), max(Counter(self.insertions))) + 1
        D = [self.deletions.count(i) for i in range(N)]
        I = [self.insertions.count(i) for i in range(N)]
        R = [i/d if d!=0 else 0 for i,d in zip(I, D)]
        fig, ax = pylab.subplots()
        ax.plot(range(N), I, marker="x", label="Insertions")
        ax.plot(range(N), D, marker="x", label="Deletions")
        ax.plot(range(N), R, "--r", label="Ratio insertions/deletions")
        ax.set_yscale("symlog")
        pylab.ylim([1, pylab.ylim()[1]])
        pylab.legend()
        pylab.grid()
        from matplotlib.ticker import MaxNLocator
        ax.xaxis.set_major_locator(MaxNLocator(integer=True))
        pylab.xlabel("Indel length", fontsize=fontsize)
        pylab.ylabel("Indel count", fontsize=fontsize)
        return I, D, R

    @_reset
    def hist_soft_clipping(self):
        """histogram of soft clipping length ignoring supplementary and
            secondary reads

        """
        from sequana import Cigar
        N = 0; M=0; C = []; F=[]
        for i,a in enumerate(self):
            c = a.cigarstring
            if c:
                C.append(Cigar(c).as_dict()['S'])
                M += 1
            else: 
                N+=1
                C.append(-1)
            F.append(a.flag)

        df = pd.DataFrame({"S":C, "F":F})
        df.query("F<32 and F!=4")['S'].hist(bins=100, log=True)
        pylab.xlabel("Soft clip length", fontsize=16)
        pylab.ylabel("#", fontsize=16)


class SAM(SAMBAMbase):
    """SAM Reader. See :class:`~samtools.bamtools.SAMBAMbase` for details"""
    def __init__(self, filename, *args):
        super(SAM, self).__init__(filename, mode="r", *args)


class CRAM(SAMBAMbase):
    """CRAM Reader. See :class:`~sequana.bamtools.SAMBAMbase` for details"""
    def __init__(self, filename, *args):
        super(CRAM, self).__init__(filename, mode="r", *args)


class BAM(SAMBAMbase):
    """BAM reader. See :class:`~sequana.bamtools.SAMBAMbase` for details"""
    def __init__(self, filename, *args):
        super(BAM, self).__init__(filename, mode="rb", *args)


class MultiBAM():
    """Convenient structure to store several BAM files

    ::

        mb = MultiBAM()
        mb.add_bam("file1.sorted.bam", group="A")
        mb.add_bam("file2.sorted.bam", group="B")

    """
    def __init__(self):
        self.bams = []
        self.tags = []
        self.groups = defaultdict(list)
        self._df = None
        self.lengths = {}

    def add_bam(self, filename, tag=None, group=None):
        from pathlib import Path
        if tag is None:
            tag = Path(filename).stem
        if group:
            self.groups[group].append(tag)
        if self.groups and group is None:
            raise ValueError('if a group was provided, it must be provided for each sample. Please provide one or create a new isntance of MultiBAM')

        bam = BAM(filename)
        self.bams.append(bam)
        self.tags.append(tag)

        # let us gather the lengths/references
        # Names must be unique and agree
        for k,v in bam.lengths.items():
            if k in self.lengths and v != self.lengths[k]:
                logger.warning("BAM reference/length incompatible with "
                    "previously loaded BAM file. Proceed but results may "
                    "be incorrect")
            self.lengths[k] = v

    def _get_df(self):
        if self._df is None:
            df = self.run()
            self._df = df
        return self._df
    df = property(_get_df)

    def run(self, exclude_secondary=True):
        from collections import Counter

        if exclude_secondary:
            data = [Counter(bam.get_df().query("flag!=4 and flag<256").rname) for bam in self.bams]
        else:
            data = [Counter(bam.get_df().query("flag!=4").rname) for bam in self.bams]

        df = pd.DataFrame(data)
        # if there is no alignments, it is equal to 0, not NA
        df = df.fillna(0)

        # let us sort the index
        df.index = self.tags
        df = df.sort_index()
        return df

    def plot_alignments_per_sample(self):
        pylab.grid(True, zorder=-1)
        if self.groups:
            N = 0
            for k, v in self.groups.items():
                data = self.df.sum(axis=1).loc[v]
                X = range(N, N+len(data))
                pylab.bar(X, data.values, label=k, ec='k', zorder=1)
                N += len(v)
            pylab.legend()
        else:
            self.df.sum(axis=1).plot(kind='bar')
        pylab.ylabel("Number of reads", fontsize=14)
        pylab.title("Number of mapped reads per sample", fontsize=14)
        pylab.tight_layout()

    def plot_alignments_per_chromosome(self):
        # on a given sample, let us keep total number of alignments (for
        # normalisation)
        S = self.df.sum(axis=1)

        # normalise and take mean number of entries per chromosomes
        mean_per_chrom = self.df.divide(S, axis=0).mean()

        # We should normalise by chromosome length. We can also simply show
        # the expected proportion for each chromosome.
        total_length = sum(list(self.lengths.values()))
        expected = [self.lengths[x]/total_length for x in mean_per_chrom.index]

        pylab.bar(range(len(expected)), expected, color='red', alpha=0.2, zorder=2)
        mean_per_chrom.plot(kind='bar', zorder=3)
        pylab.grid(True, zorder=1)
        pylab.legend(['Expected', 'Observed'])
        pylab.ylabel('Proportion alignements per chromosome', fontsize=14)
        return mean_per_chrom, expected

class Alignment(object):
    """Helper class to retrieve info about Alignment

    Takes an alignment as read by :class:`BAM` and provides a simplified version
    of pysam.Alignment class.

    ::

        >>> from sequana.bamtools import Alignment
        >>> from sequana import BAM, sequana_data
        >>> b = BAM(sequana_data("test.bam"))
        >>> segment = next(b)
        >>> align = Alignment(segment)
        >>> align.as_dict()
        >>> align.FLAG
        353

    The original data is stored in hidden attribute :attr:`_data` and the
    following values are available as attributes or dictionary:


    * QNAME: a query template name. Reads/segment having same QNAME come from the
      same template. A QNAME set to `*` indicates the information is unavailable.
      In a sam file, a read may occupy multiple alignment
    * FLAG: combination of bitwise flags. See :class:`SAMFlags`
    * RNAME: reference sequence
    * POS
    * MAPQ: mapping quality if segment is mapped. equals -10 log10 Pr
    * CIGAR: See :class:`sequana.cigar.Cigar`
    * RNEXT: reference sequence name of the primary alignment of the NEXT read
      in the template
    * PNEXT: position of primary alignment
    * TLEN: signed observed template length
    * SEQ: segment sequence
    * QUAL: ascii of base quality


    """
    def __init__(self, alignment):
        """.. rubric:: constructor

        :param alignment: alignment instance from :class:`BAM`


        """
        self._data = alignment
        d = self.as_dict()
        for key in d.keys():
            setattr(self, key, d[key])

    def as_dict(self):
        d = {}
        s = self._data
        d['QNAME'] = s.qname
        d['FLAG'] = s.flag
        d['RNAME'] = s.rname
        d['POS'] = s.pos
        d['MAPQ'] = s.mapq
        d['CIGAR'] = s.cigar
        d['PNEXT'] = s.pnext
        d['RNEXT'] = s.rnext
        d['TLEN'] = s.tlen
        d['SEQ'] = s.seq
        d['QUAL'] = s.qual
        return d


class SAMFlags(object):
    """Utility to extract bits from a SAM flag

    .. doctest::

        >>> from sequana import SAMFlags
        >>> sf = SAMFlags(257)
        >>> sf.get_flags()
        [1, 256]


    You can also print the bits and their description::

        print(sf)

    ======= ====================================================================
    bit     Meaning/description
    ======= ====================================================================
    0       mapped segment
    1       template having multiple segments in sequencing
    2       each segment properly aligned according to the aligner
    4       segment unmapped
    8       next segment in the template unmapped
    16      SEQ being reverse complemented
    32      SEQ of the next segment in the template being reverse complemented
    64      the first segment in the template
    128     the last segment in the template
    256     secondary alignment
    512     not passing filters, such as platform/vendor quality controls
    1024    PCR or optical duplicate
    2048    supplementary alignment
    ======= ====================================================================

    :reference: http://samtools.github.io/hts-specs/SAMv1.pdf
    """
    def __init__(self, value=4095):
        self.value = value
        self._flags = {
            0: "segment mapped",
            1: "template having multiple segments in sequencing",
            2: "each segment properly aligned according to the aligner",
            4: "segment unmapped",
            8: "next segment in the template unmapped",
            16: "SEQ being reverse complemented",
            32: "SEQ of the next segment in the template being reverse complemented",
            64: "the first segment in the template",
            128: "the last segment in the template",
            256: "secondary alignment",
            512: "not passing filters, such as platform/vendor quality controls",
            1024: "PCR or optical duplicate",
            2048: "supplementary alignment"}

    def get_meaning(self):
        """Return all description sorted by bit """
        return [self._flags[k] for k in sorted(self._flags.keys())]

    def get_flags(self):
        """Return the individual bits included in the flag"""
        flags = []
        for this in sorted(self._flags.keys()):
            if self.value & this:
                flags.append(this)
        return flags

    def __str__(self):
        txt = ""
        for this in sorted(self._flags.keys()):
            if self.value & this:
                txt += "%s: %s\n" % (this, self._flags[this])
        return txt


class CS(dict):
    """Interpret CS tag from SAM/BAM file tag

    ::

        >>> from sequana import CS
        >>> CS('-a:6-g:14+g:2+c:9*ac:10-a:13-a')
        {'D': 3, 'I': 2, 'M': 54, 'S': 1}

    When using some mapper, CIGAR are stored in another format called CS, which
    also includes the substitutions. See minimap2 documentation for details.
    """
    def __init__(self, tag):
        self.tag = tag
        d = self._scan()
        for k,v in d.items():
            self[k] = v

    def _scan(self):
        d = {"M":0, "I":0, "D":0, "S":0}
        current = ":"  # this is just to start the loop with a key (set to 0)
        number = "0"

        for c in self.tag:
            if c in ":+-*":
                if current == ":":
                    d["M"] += int(number)
                elif current == "+":
                    d["I"] += len(number)
                elif current == "-":
                    d["D"] += len(number)
                elif current == "*":
                    d["S"] += len(number)
                current = c
                number = ""
            else: # a letter or number
                number += c
        # last one
        if current == ":":
            d["M"] += int(number)
        elif current == "+":
            d["I"] += len(number)
        elif current == "-":
            d["D"] += len(number)
        elif current == "*":
            d["S"] += len(number)

        assert d['S'] % 2 == 0
        d['S'] //= 2
        return d


