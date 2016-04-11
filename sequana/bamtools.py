"""Tools to manipulate BAM/SAM files

.. autosummary::

    Alignment
    BAM
    SAMFlags

.. note:: BAM being the compressed version of SAM files, we do not
    implement any function to read SAM files. We strongly encourage
    developers to convert their SAM to BAM

"""
import os
import pandas as pd
import pylab
import pysam

"""
#http://www.acgt.me/blog/2014/12/16/understanding-mapq-scores-in-sam-files-does-37-42#
#http://biofinysics.blogspot.fr/2014/05/how-does-bowtie2-assign-mapq-scores.html
#https://gitlab.univ-nantes.fr/a-slide/ContaVect/blob/9a411abfa720064c205c5f6c811afdfea206ed12/pyDNA/pySamTools/Bam.py
# pysamtools
# pysam uses htlib behing the scene and is very fast
# pysam works great for BAM file but with SAM, it needs to read the file after
# each compete iteration, which is not very useful

Interesting commands::

    samtools flagstat contaminant.bam
    samtools stats contaminant.bam
"""


# simple decorator ro rewind the BAM file
from functools import wraps
def seek(f):
    @wraps(f)
    def wrapper(*args, **kargs):
        #args[0] is the self of the method
        # For BAM files only
        args[0].reset()
        return f(*args, **kargs)
    return wrapper


class BAM(pysam.AlignmentFile):
    """BAM data structure

    This is built on top of pysam package and provides functions to retrieve
    useful information or plot various statistical analysis.

    This BAM class can also read a SAM file but some functionalities will not
    work. Besides, you need to create new instances each time a method is
    called. This is inherent to pysam implementation. Python2.7 and 3.5 also
    behave differently and we would recommend the py35 version. For instance,
    :meth:`to_fastq` would work only for that version;

    .. doctest::

        >>> from sequana import BAM, sequana_data
        >>> b = BAM(sequana_data("test.bam", "testing"))
        >>> len(b)
        1000


    .. note:: Once you loop over this data structure,  you must call
        :meth`reset` to rewind the loop. The methods implemented in this data
        structure takes care of that for you thanks to a decorator called seek.
        So if you want to call an iterator yourself you must use the
        :meth:`reset` method yourself. In particular, if you want to use the
        next() funtion.



    """
    def __init__(self, filename,  mode='rb', *args):
        # The mode rb means read-only (r) and that the format is BAM or SAM (b)

        # super()  works for py27 but not py35 probably a missing __init__  or __new__ in
        # AlignmentFile class. See e.g., stackoverflow/questions/
        # 26653401/typeerror-object-takes-no-parameters-but-only-in-python-3

        #self._filename = filename
        pysam.AlignmentFile.__init__(filename, mode=mode, *args)

        self._filename = filename
        # the BAM format can be rewinded but not SAM. This is a pain since one
        # need to reload the instance after each operation. With the BAM, we can
        # do that automatically. We therefore enfore the BAM format.

        # Another issue is that in PY3, there is an attribute **format** that
        # tells us if the input is BAM or SAM but this does not work in PY2
        # where the attribute is missing...

        # another know issue is that we cannot use the method is_bam trustfully.
        # Using AlignmentFile one must provide the proper mode 'e.g. rb will set
        # is_bam to True while 'r' only will set it to False. However, it seems
        # there is not sanity check inside the file. Besides, using
        # AlignmentFile.__init__ instead of the class itself seems to haev also
        # some side effects.

        # Let us save the length (handy and use in other places).
        # If is is a SAM file, the rewind does not work and calling it again wil
        # return 0. This may give us a hint that it is a SAM file
        self.N = len(self)

        # running a second time the len() should return the correct answer with
        # BAM files but the SAM will not work and return 0
        if len(self) == 0:
            raise ValueError("Convert your SAM file to a BAM file please")

        # since we called len(), let us reset the iterator so that if the user
        # decide to use the BAM data structure, it will work at least once as
        # expected by the pysam API
        self.reset()

    @seek
    def get_read_names(self):
        """Return the reads' names"""
        names = [this.qname for this in self]
        return names

    @seek
    def iter_unmapped_reads(self):
        """Return an iterator on the reads that are unmapped"""
        unmapped = (this.qname for this in self if this.is_unmapped)
        return unmapped

    @seek
    def iter_mapped_reads(self):
        """Return an iterator on the reads that are mapped"""
        mapped = (this.qname for this in self if this.is_unmapped is False)
        return mapped

    @seek
    def __len__(self):
        N = len([x for x in self])
        return N

    @seek
    def get_flags(self):
        """Return flags of all reads as a list


        .. seealso:: :meth:`get_flags_as_df`, :class:`SAMFlags`"""
        flags = [s.flag for s in self]
        return flags

    def get_stats(self):
        """Return basic stats about the reads

        :return: dictionary with basic stats:

            - total_reads : number reads
            - mapped_reads : number of mapped reads
            - unmapped_reads : number of unmapped
            - contamination [%]: mapped / unmapped

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
        d['total_reads'] = len(list(self.iter_unmapped_reads()))
        d['mapped_reads'] = len(list(self.iter_mapped_reads()))
        d['unmapped_reads'] = len(list(self.iter_unmapped_reads()))
        d['contamination [%]'] = float(d['mapped_reads']) /float(d['unmapped_reads'])
        d['contamination [%]'] *= 100
        return d

    def get_full_stats_as_df(self):
        """Return a dictionary with full stats about the BAM file

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

    def get_flags_as_df(self):
        """Returns flags as a dataframe

        .. doctest::

            >>> from sequana import BAM, sequana_data
            >>> b = BAM(sequana_data('test.bam', "testing"))
            >>> df = b.get_flags_as_df()
            >>> df.sum()
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
        flags = self.get_flags()
        data = [(this, [flag&this for flag in flags])
            for this in [1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048]]
        import pandas as pd
        df = pd.DataFrame(dict(data))
        df = df > 0
        return df

    def plot_bar_flags(self, logy=True, fontsize=16, filename=None):
        """

        .. plot::
            :include-source:

            from sequana import BAM, sequana_data
            b = BAM(sequana_data('test.bam', "testing"))
            b.plot_bar_flags()


        """
        df = self.get_flags_as_df()
        df = df.sum()
        pylab.clf()
        if logy is True:
            df.plot(kind='bar', logy=logy, grid=True)
        else:
            df.plot(kind='bar', grid=True)
        pylab.xlabel("flags", fontsize=fontsize)
        pylab.ylabel("count", fontsize=fontsize)
        pylab.tight_layout()
        if filename:
            pylab.savefig(filename)

    @seek
    def to_fastq(self, filename):
        """Export the BAM to FastQ format

        .. warning:: to be tested
        .. todo:: comments from original reads ?
        """
        with open(filename, "w") as fh:
            for i, this in enumerate(self):
                read = this.qname + "\t" + "???\n"
                read += this.seq + "\n"
                read += "+\n"
                read += this.qual + "\n"
                #if i != self.N-1:
                #    read += "\n"
                fh.write(read)

    @seek
    def get_mapq_as_df(self):
        """Return dataframe with mapq for each read"""
        df = pd.DataFrame({'mapq': [this.mapq for this in self]})
        return df

    def plot_bar_mapq(self, fontsize=16, filename=None):
        """Plots bar plots of the MAPQ (quality) of alignments

            .. plot::
                :include-source:

                from sequana import BAM, sequana_data
                b = BAM(sequana_data('test.bam', "testing"))
                b.plot_bar_mapq()

        """
        df = self.get_mapq_as_df()
        df.plot(kind='hist', bins=60, legend=False, grid=True, logy=True)
        pylab.xlabel("MAPQ", fontsize=fontsize)
        try:
            # This may raise issue on MAC platforms
            pylab.tight_layout()
        except:
            pass
        if filename:
            pylab.savefig(filename)


class Alignment(object):
    """Helper class to retrieve info about Alignment

    Takes an alignment as read by :class:`BAM` and provide simplified version
    of pysam.Alignment class.

    ::

        >>> from sequana.bamtools import Alignment
        >>> from sequana import BAM, sequana_data
        >>> b = BAM(sequana_data("test.bam", "testing"))
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
    * CIGAR:
    * RNEXT: reference sequence name of the primary alignment of the NEXT read
      in the template
    * PNEXT: position of primary alignment
    * TLEN: signed observed template length
    * SEQ: segment sequence
    * QUAL: ascii of base quality

    .. note:: A segment is a contiguous sequence. A read is a sequence that
        may consist of multiple segments.


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
    2048    supplementary alignme
    ======= ====================================================================

    """
    def __init__(self, value=4095):
        self.value = value
        self._flags = {
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
            2048: "supplementary alignme"}

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


