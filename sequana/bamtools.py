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
    SAM
    SAMFlags

.. note:: BAM being the compressed version of SAM files, we do not
    implement any functionalities related to SAM files. We strongly encourage
    developers to convert their SAM to BAM.

"""
import os
import json
from collections import Counter
from collections import OrderedDict

from sequana.lazy import pandas as pd
from sequana.lazy import numpy as np
from sequana.lazy import pylab

import pysam
from sequana import jsontool, logger
logger.name = __name__
"""
#http://www.acgt.me/blog/2014/12/16/understanding-mapq-scores-in-sam-files-does-37-42#
#http://biofinysics.blogspot.fr/2014/05/how-does-bowtie2-assign-mapq-scores.html
#https://gitlab.univ-nantes.fr/a-slide/ContaVect/blob/9a411abfa720064c205c5f6c811afdfea206ed12/pyDNA/pySamTools/Bam.py

Interesting commands::

    samtools flagstat contaminant.bam
"""

__all__ = ['BAM','Alignment', 'SAMFlags', "CS", "SAM", "CRAM"]


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
    def get_df(self, max_align=-1):
        flags = []
        starts = []
        ends = []
        mapqs = []
        refnames = []
        querynames = []
        querylengths = []
        for i, a in enumerate(self._data):
            flags.append(a.flag)
            starts.append(a.reference_start)
            ends.append(a.reference_end)
            mapqs.append(a.mapq)
            try:refnames.append(a.reference_name)
            except:refnames.append(-1)
            querynames.append(a.query_name)
            querylengths.append(a.query_length)
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
        I, D, M, L, mapq, flags, NM = [], [], [], [], [], [], []
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

        try:
            S = np.array(S)
            C = 1 - (I + D + S)/(S + I + D + M)
            logger.info("computed Concordance based on minimap2 --cs option")
        except:
            logger.info("computed Concordance based on standard CIGAR information using INDEL and NM tag")
            computed_S = NM - D - I
            C = 1 - (I + D + computed_S)/(computed_S + I + D + M)

        df = pd.DataFrame([C, L, I, D, M, mapq, flags, NM, S])
        df = df.T
        df.columns = ["concordance", 'length', "I", "D", "M", "mapq", "flags", "NM", "mismatch"]
        return df

    def __iter__(self):
        return self

    def __next__(self):
        return next(self._data)

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
            if this.pos < pos:
                self._sorted = False
                return False
            pos = this.pos
        self._sorted = True
        return self._sorted
    is_sorted = property(_get_is_sorted, doc="return True if the BAM is sorted")

    def get_full_stats_as_df(self):
        """Return a dictionary with full stats about the BAM/SAM file

        The index of the dataframe contains the flags. The column contains
        the counts.

        ::

            >>> from sequana import BAM, sequana_data
            >>> b = BAM(sequana_data("test.bam"))
            >>> df = b.get_full_stats_as_df()
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
        for read in self:
            self._count_item(mapq_dict, read.mapq)
            self._count_item(flag_dict, read.flag)
            if read.is_unmapped is False:
                self._count_item(read_length_dict, read.reference_length)
            try:mean_qualities.append(pylab.mean(read.query_qualities))
            except:mean_qualities.append(read.query_qualities)
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

        .. warning:: works only for BAM files. Use :meth:`get_full_stats_as_df`
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
    def get_mapq_as_df(self):
        """Return dataframe with mapq for each read"""
        df = pd.DataFrame({'mapq': [this.mapq for this in self]})
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
        data = [(f.seq.count("C") + f.seq.count('G')) / len(f.seq)*100. for f in self]
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
        df = pd.DataFrame(qualities)
        from sequana.viz.boxplot import Boxplot
        bx = Boxplot(df)
        try:
            bx.plot(ax=ax)
        except:
            bx.plot()

    def _set_alignments(self):
        # this scans the alignments once for all
        self.alignments = [this for this in self]

    @_reset
    def _set_coverage(self):
        try: self.alignments
        except: self._set_alignments()

        reference_start = [this.reference_start for this in self.alignments]
        reference_end = [this.reference_end for this in self.alignments]
        N = max([this for this in reference_end if this])

        self.coverage = np.zeros(N)
        for x, y in zip(reference_start, reference_end):
            if y and x>=0 and y>=0: self.coverage[x:y] += 1
            else: pass

    @_reset
    def _set_indels(self):
        try: self.alignments
        except: self._set_alignments()

        self.insertions = []
        self.deletions = []
        for this in self.alignments:
            if this.cigarstring:
                if "I" in this.cigarstring:
                    self.insertions.extend([x[1] for x in this.cigartuples if x[0] == 1])
                if "D" in this.cigarstring:
                    self.deletions.extend([x[1] for x in this.cigartuples if x[0] == 2])

    def plot_coverage(self):
        """Please use :class:`GenomeCov` for more sophisticated
        tools to plot the genome coverage

        .. plot::
            :include-source:

            from sequana import sequana_data, BAM
            b = BAM(sequana_data("measles.fa.sorted.bam"))
            b.plot_coverage()

        """
        try: self.coverage
        except: self._set_coverage()
        pylab.plot(self.coverage)
        pylab.xlabel("Coverage")

    def hist_coverage(self, bins=100):
        """

        .. plot::
            :include-source:

            from sequana import sequana_data, BAM
            b = BAM(sequana_data("measles.fa.sorted.bam"))
            b.hist_coverage()
        """
        try: self.coverage
        except: self._set_coverage()
        pylab.hist(self.coverage, bins=bins)
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
    """SAM Reader. See :class:`~samtools.bamtools.SAMBAMBase` for details"""
    def __init__(self, filename, *args):
        super(SAM, self).__init__(filename, mode="r", *args)

class CRAM(SAMBAMbase):
    """CRAM Reader. See :class:`~sequana.bamtools.SAMBAMBase` for details"""
    def __init__(self, filename, *args):
        super(CRAM, self).__init__(filename, mode="r", *args)


class BAM(SAMBAMbase):
    """BAM reader. See :class:`~sequana.bamtools.SAMBAMBase` for details"""
    def __init__(self, filename, *args):
        super(BAM, self).__init__(filename, mode="rb", *args)





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


