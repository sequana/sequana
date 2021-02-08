# -*- coding: utf-8 -*-
#
#  This file is part of Sequana software
#
#  Copyright (c) 2016-2017 - Sequana Development Team
#
#  File author(s):
#      Thomas Cokelaer <thomas.cokelaer@pasteur.fr>
#      Mélissa Cardon <melissa.cardon@pasteur.fr>,
#
#  Distributed under the terms of the 3-clause BSD license.
#  The full license is in the LICENSE file, distributed with this software.
#
#  website: https://github.com/sequana/sequana
#  documentation: http://sequana.readthedocs.io
#
##############################################################################
"""Pacbio QC and stats"""
import os
import collections
import json
import random

from sequana.lazy import pylab
from sequana.lazy import numpy as np
from sequana.lazy import pandas as pd
import pysam

import colorlog
logger = colorlog.getLogger(__name__)


from sequana.summary import Summary


__all__ = ["PacbioMappedBAM", "PacbioSubreads", "PBSim", "BAMSimul", "Barcoding"]


# from pbcore.io import openAlignmentFile
# b = openAlignmentFile(filename)
# len(b) is instantaneous IF a bai is created using
# samtools index -b filename.bam filename.bai


class HistCumSum(object):

    def __init__(self, data, grid=True, fontsize=16, xlabel="", ylabel="",
                 title=""):
        self.data = data
        self.fontsize = fontsize
        self.xlabel = xlabel
        self.ylabel = ylabel
        self.title = title
        self.grid = grid

    def plot(self, bins=80, rwidth=0.8, **kwargs):
        pylab.clf()
        Y, X, _ = pylab.hist(self.data, bins=bins, rwidth=rwidth, **kwargs)

        pylab.xlabel(self.xlabel, fontsize=self.fontsize)
        pylab.ylabel(self.ylabel, fontsize=self.fontsize)

        """self.Y = Y
        self.X = X

        ax_twin = pylab.gca().twinx()

        shift = (X[1] - X[0]) / 2

        ax_twin.plot(X[0:-1]- shift, len(self.data) - pylab.cumsum(Y), "k")
        ax_twin.set_ylim(bottom=0)
        pylab.ylabel("CDF", fontsize=self.fontsize)
        """
        pylab.grid(self.grid)
        pylab.title(self.title)
        try: pylab.tight_layout()
        except:pass


class PacbioBAMBase(object):
    """Base class for Pacbio BAM files


    """
    def __init__(self, filename):
        """

        :param str filename: input BAM file

        """
        self.filename = filename
        self.data = pysam.AlignmentFile(filename, check_sq=False)
        self._df = None
        self._nb_pass = None
        self.sample_name = os.path.basename(self.filename)

    def __len__(self):
        return len(self.df)

    def __str__(self):
        return "Length: {}".format(len(self))

    def reset(self):
        self.data.close()
        self.data = pysam.AlignmentFile(self.filename, check_sq=False)

    def _to_fastX(self, mode, output_filename, threads=2):
        """

        :param mode: fastq or fasta

        """
        # for now, we use samtools
        # can use bamtools as well but as long and output 10% larger (sequences
        # are split on 80-characters length)
        from snakemake import shell
        cmd = "samtools %s  -@ %s %s > %s" % (mode, threads,
            self.filename, output_filename)
        logger.info("Please be patient")
        logger.info("This may be long depending on your input data file: ")
        logger.info("typically, a minute per  500,000 reads")
        shell(cmd)
        logger.info("done")

    def to_fastq(self, output_filename, threads=2):
        """Export BAM reads into FastQ file"""
        self._to_fastX("fastq", output_filename, threads=threads)

    def to_fasta(self, output_filename, threads=2):
        """Export BAM reads into a Fasta file

        :param output_filename: name of the output file (use .fasta extension)
        :param int threads: number of threads to use

        .. note:: this executes a shell command based on samtools

        .. warning:: this takes a few minutes for 500,000 reads

        """
        self._to_fastX("fasta", output_filename, threads=threads)

    def hist_GC(self, bins=50, alpha=0.5, hold=False, fontsize=12,
                grid=True, xlabel="GC %", ylabel="#", label="",title=None):
        """Plot histogram GC content

        :param int bins: binning for the histogram
        :param float alpha: transparency of the histograms
        :param bool hold:
        :param int fontsize: fontsize of the x and y labels and title.
        :param bool grid: add grid or not
        :param str xlabel:
        :param str ylabel:
        :param str label: label of the histogram (for the legend)
        :param str title:

        .. plot::
            :include-source:

            from sequana.pacbio import PacbioSubreads
            from sequana import sequana_data
            b = PacbioSubreads(sequana_data("test_pacbio_subreads.bam"))
            b.hist_GC()

        """
        mean_GC =  np.mean(self.df.loc[:,'GC_content'])

        # set title if needed
        if title is None:
            title = "GC %%  \n Mean GC : %.2f" %(mean_GC)

        # histogram GC percent
        if hold is False:
            pylab.clf()
        pylab.hist(self.df.loc[:,'GC_content'], bins=bins,
            alpha=alpha, label=label + ", mean : " + str(round(mean_GC, 2))
            + ", N : " + str(len(self)))
        pylab.xlabel(xlabel, fontsize=fontsize)
        pylab.ylabel(ylabel, fontsize=fontsize)
        pylab.title(title, fontsize=fontsize)
        if grid is True:
            pylab.grid(True)
        pylab.xlim([0, 100])
        try: pylab.tight_layout()
        except:pass

    def plot_GC_read_len(self, hold=False, fontsize=12, bins=[200, 60],
                grid=True, xlabel="GC %", ylabel="#", cmap="BrBG"):
        """Plot GC content versus read length

        :param bool hold:
        :param int fontsize: for x and y labels and title
        :param bins: a integer or tuple of 2 integers to specify
            the binning of the x and y 2D histogram.
        :param bool grid:
        :param str xlabel:
        :param str ylabel:

        .. plot::
            :include-source:

            from sequana.pacbio import PacbioSubreads
            from sequana import sequana_data
            b = PacbioSubreads(sequana_data("test_pacbio_subreads.bam"))
            b.plot_GC_read_len(bins=[10, 10])

        """
        from sequana.viz import Hist2D
        mean_len =  np.mean(self.df.loc[:,'read_length'])
        mean_GC =  np.mean(self.df.loc[:,'GC_content'])

        if hold is False:
            pylab.clf()


        data = self.df.loc[:,['read_length','GC_content']].dropna()
        h = Hist2D(data)
        res = h.plot(bins=bins, contour=False, norm='log', Nlevels=6, cmap=cmap)
        pylab.xlabel("Read length", fontsize=fontsize)
        pylab.ylabel("GC %", fontsize=fontsize)
        pylab.title("GC %% vs length \n Mean length : %.2f , Mean GC : %.2f" %
            (mean_len, mean_GC), fontsize=fontsize)
        pylab.ylim([0, 100])
        if grid is True:
            pylab.grid(True)

    def hist_read_length(self, bins=80, alpha=0.5, hold=False, fontsize=12,
                grid=True, xlabel="Read Length", ylabel="#", label="",
                title=None, logy=False,  ec="k", hist_kwargs={}):
        """Plot histogram Read length

        :param int bins: binning for the histogram
        :param float alpha: transparency of the histograms
        :param bool hold:
        :param int fontsize:
        :param bool grid:
        :param str xlabel:
        :param str ylabel:
        :param str label: label of the histogram (for the legend)
        :param str title:

        .. plot::
            :include-source:

            from sequana.pacbio import PacbioSubreads
            from sequana import sequana_data
            b = PacbioSubreads(sequana_data("test_pacbio_subreads.bam"))
            b.hist_read_length()

        """
        mean_len =  np.mean(self.df.loc[:,'read_length'])

        # set title if not provided
        if title is None:
            title = "Read length  \n Mean length : %.2f" %(mean_len)

        if hold is False:
            pylab.clf()

        hist = HistCumSum(self.df.loc[:,'read_length'], fontsize=fontsize,
                    grid=grid)
        hist.title = title
        hist.xlabel = xlabel
        hist.ylabel = ylabel
        hist.plot(bins=bins, alpha=alpha, edgecolor=ec,
            label=  "%s, mean : %.0f, N : %d" % (label, mean_len, len(self)),
            log=logy, **hist_kwargs)
        pylab.gca().set_ylim(bottom=0)
        pylab.gca().set_xlim(left=0)


class PacbioSubreads(PacbioBAMBase):
    """BAM reader for Pacbio (reads)

    You can read a file as follows::

        from sequana.pacbio import Pacbiosubreads
        from sequana import sequana_data
        filename = sequana_data("test_pacbio_subreads.bam")
        b = PacbioSubreads(filename)

    A summary of the data is stored in the attribute :attr:`df`. It contains
    information such as the length of the reads, the ACGT content, the GC content.

    Several plotting methods are available. For instance, :meth:`hist_snr`.


    The BAM file used to store the Pacbio reads follows the BAM/SAM
    specification. Note that the sequence read are termed query, a subsequence
    of an entire Pacbio ZMW read ( a subread), which is basecalls from a single
    pass of the insert DNA molecule.

    In general, only a subsequence of the query will align to the
    reference genome, and that subsequence is referred to as the aligned query.

    When introspecting the aligned BAM file, the extent of the query in ZMW read
    is denoted as [qStart, qEnd) and the extent of the aligned subinterval
    as [aStart, aEnd). The following graphic illustrates these intervals::

              qStart                         qEnd
        0         |  aStart                aEnd  |
        [--...----*--*---------------------*-----*-----...------)  < "ZMW read" coord. system
                  ~~~----------------------~~~~~~                  <  query; "-" =aligning subseq.
        [--...-------*---------...---------*-----------...------)  < "ref." / "target" coord. system
        0            tStart                tEnd

    In the BAM files, the qStart, qEnd are contained in the qs and qe tags, (and
    reflected in the QNAME); the bounds of the aligned query in the ZMW read can be
    determined by adjusting qs and qe by the number of soft-clipped bases at the
    ends of the alignment (as found in the CIGAR).

    See also the comments in the code for other tags.

    :reference: http://pacbiofileformats.readthedocs.io/en/3.0/BAM.html


    """
    def __init__(self, filename, sample=0):
        """.. rubric:: Constructor

        :param str filename: filename of the input pacbio BAM file. The content
            of the BAM file is not the ouput of a mapper. Instead, it is the
            output of a Pacbio (Sequel) sequencing (e.g., subreads).
        :param int sample: for sample, you can set the number of subreads to
            read (0 means read all subreads)
        """
        super(PacbioSubreads, self).__init__(filename)
        self._sample = sample

    def _get_df(self):
        # When scanning the BAM, we can extract the length, SNR of ACGT (still
        # need to know how to use it). The GC content (note there is no
        # ambiguity so no S character). The ZMW. Also, from the tags we could
        # get more

        # In each alignement, there are lots of information to retrieve.
        # One could for instance introspect the tags.
        # - cx: subread local context flags
        # - ip: vector of length qlen from 0 to 250. This is the IPD (raw frames
        # or codec V1)
        # - np: number of passes (1 for subread, variable for CCS)
        # - pw: vector of length qlen from 0 to 128? This is the PulseWidth (raw
        # frames or codec V1)
        # - qs: 0-based start of query in the ZMW read (absent in CCS)
        # - qe: 0-based end of query in the ZMW read (absent in CCS)
        # - zm: position/ID of the ZMW
        # - sn: list of ACGT SNRs. A, C, G, T in that order
        # - rq: float encoding exepted accuracy

        # - dq: DeletionQV
        # - dt: deletion Tag
        # - iq: insertionQV
        # - mq: mergeQV
        # - sq: substituionQV
        # - st: substituion tag
        # - RG: ?

        # See http://pacbiofileformats.readthedocs.io/en/3.0/BAM.html
        if self._df is None:
            logger.info("Scanning input file. Please wait")
            self.reset()
            N = 0

            all_results = []
            # This takes 60%  of the time...could use cython ?
            for i, read in enumerate(self.data):
                tags = dict(read.tags)
                res = []
                # count reads
                N += 1
                if (N % 10000) == 0:
                    logger.info("Read %d sequences" %N)

                # res[0] = read length
                res.append(read.query_length) # also stored in tags["qe"] - tags["qs"]
                res.append(read.reference_length) # also stored in tags["qe"] - tags["qs"]

                # collections.counter is slow, let us do it ourself
                if read.query_length and read.query_sequence:
                    res.append( 100. / read.query_length * sum(
                        [read.query_sequence.count(letter) for letter in "CGcgSs"]))
                else:
                    res.append(None)

                # res[1:4] contains SNR  stored in tags['sn'] in the order A, C, G, T
                try:
                    snr = list(tags['sn'])
                except:
                    snr = [None] * 4
                res = res + snr

                # res[6] = ZMW name, also stored in tags["zm"]
                try:
                    res.append(int(read.qname.split('/')[1]))
                except: # simulated data may not have the ZMW info, in which
                        #case, we store just a unique ID
                    res.append(i)

                # aggregate results
                all_results.append(res)

                if self._sample and N >= self._sample:
                    break

            self._df = pd.DataFrame(all_results,
                columns=['read_length', "reference_length", 'GC_content',
                            'snr_A','snr_C','snr_G','snr_T','ZMW'])

            # populate the nb passes from the ZMW
            grouped = self._df.groupby("ZMW")
            agg = grouped.agg({"read_length": len})

            ZMW = self._df.ZMW.unique()
            aa = list(pylab.flatten([[agg.loc[this][0]] * agg.loc[this][0] for this in ZMW]))
            self._df['nb_passes'] = aa
            self._df['nb_passes'] -= 1 # nb passes starts at 0

            self.reset()
        return self._df
    df = property(_get_df)

    def _get_stats(self):
        # cast to allows json dump
        data =  self.df.read_length.describe().to_dict()
        data['nb_bases'] = int(self.df.read_length.sum())
        data['nb_reads'] = len(self.df)
        data['mean_GC'] = float(self.df.GC_content.mean())
        data['CCS_nb_reads'] = self.get_number_of_ccs()
        data['CCS_mean_passes'] = self.get_number_of_ccs()
        return data
    stats = property(_get_stats, doc="return basic stats about the read length")

    def stride(self, output_filename, stride=10, shift=0, random=False):
        """Write a subset of reads to BAM output

        :param str output_filename: name of output file
        :param int stride: optionnal, number of reads to read to output one read
        :param int shift: number of reads to ignore at the begining of input file
        :param bool random: if True, at each step the read to output is randomly selected
        """
        assert output_filename != self.filename, \
            "output filename should be different from the input filename"
        self.reset()
        with pysam.AlignmentFile(output_filename,"wb", template=self.data) as fh:
            if random:
                shift = np.random.randint(stride)

            for i, read in enumerate(self.data):
                if (i + shift) % stride == 0:
                    fh.write(read)
                    if random:
                        shift = np.random.randint(stride)

    def random_selection(self, output_filename, nreads=None,
            expected_coverage=None, reference_length=None, read_lengths=None):
        """Select random reads

        :param nreads: number of reads to select randomly. Must be less than
            number of available reads in the orignal file.
        :param expected_coverage:
        :param reference_length:

        if expected_coverage and reference_length provided, nreads is replaced
        automatically.

        .. note:: to speed up computation (if you need to call random_selection
            many times), you can provide the mean read length manually
        """
        assert output_filename != self.filename, \
            "output filename should be different from the input filename"

        if read_lengths is None:
            self.reset()
            read_lengths = [read.query_length for i, read in enumerate(self.data)]

        N = len(read_lengths)

        if expected_coverage and reference_length:
            mu = pylab.mean(read_lengths)
            nreads = int(expected_coverage * reference_length / mu)

        assert nreads < N, "nreads parameter larger than actual Number of reads"
        selector = random.sample(range(N), nreads)
        logger.info("Creating a pacbio BAM file with {} reads".format(nreads))

        with pysam.AlignmentFile(output_filename,"wb", template=self.data) as fh:
            self.reset()
            for i, read in enumerate(self.data):
                if i in selector:
                    fh.write(read)

    def summary(self):
        summary = {"name": "sequana_summary_pacbio_qc"}
        summary["read_stats"] = self.stats.copy()
        summary["mean_gc"] = float(np.mean(self._df.loc[:,'GC_content']))
        a, b = np.histogram(self._df.loc[:,'GC_content'], bins=100)
        summary['hist_gc'] = {"Y": a.tolist(), "X": b.tolist()}
        a, b =  np.histogram(self._df['read_length'],100)
        summary['hist_read_length'] = {"Y": a.tolist(), "X": b.tolist()}
        return summary

    def save_summary(self, filename):
        summary = self.summary()
        with open(filename, "w") as fh:
            json.dump(summary, fh, indent=4, sort_keys=True)

    def filter_length(self, output_filename, threshold_min=0,
        threshold_max=np.inf):
        """Select and Write reads within a given range

        :param str output_filename: name of output file
        :param int threshold_min: minimum length of the reads to keep
        :param int threshold_max: maximum length of the reads to keep

        """
        assert threshold_min < threshold_max
        assert output_filename != self.filename, \
            "output filename should be different from the input filename"
        self.reset()
        with pysam.AlignmentFile(output_filename,  "wb", template=self.data) as fh:
            for read in self.data:
                if ((read.query_length > threshold_min) & (read.query_length < threshold_max)):
                    fh.write(read)

    def hist_snr(self, bins=50, alpha=0.5, hold=False, fontsize=12,
                grid=True, xlabel="SNR", ylabel="#",title="", clip_upper_SNR=30):
        """Plot histogram of the ACGT SNRs for all reads

        :param int bins: binning for the histogram. Note that the range starts
            at 0 and ends at clip_upper_SNR
        :param float alpha: transparency of the histograms
        :param bool hold:
        :param int fontsize:
        :param bool grid:
        :param str xlabel:
        :param str ylabel:
        :param str title:

        .. plot::
            :include-source:

            from sequana.pacbio import PacbioSubreads
            from sequana import sequana_data
            b = PacbioSubreads(sequana_data("test_pacbio_subreads.bam"))
            b.hist_snr()

        """
        if self._df is None:
            self._get_df()

        # old pacbio format has no SNR stored
        if len(self._df['snr_A'].dropna()) == 0:
            # nothing to plot
            from sequana import sequana_data
            pylab.clf()
            pylab.imshow(pylab.imread(sequana_data("no_data.jpg")))
            pylab.gca().axis('off')
            return

        if hold is False:
            pylab.clf()

        maxSNR = 0
        for letter in "ACGT":
            m = self._df.loc[:,"snr_{}".format(letter)].max()
            if m > maxSNR:
                maxSNR = m

        if maxSNR > clip_upper_SNR:
            maxSNR = clip_upper_SNR

        bins = pylab.linspace(0, maxSNR, bins)

        pylab.hist(self._df.loc[:,'snr_A'].clip(upper=maxSNR), alpha=alpha, label="A", bins=bins)
        pylab.hist(self._df.loc[:,'snr_C'].clip(upper=maxSNR), alpha=alpha, label="C", bins=bins)
        pylab.hist(self._df.loc[:,'snr_G'].clip(upper=maxSNR), alpha=alpha, label="G", bins=bins)
        pylab.hist(self._df.loc[:,'snr_T'].clip(upper=maxSNR), alpha=alpha, label="T", bins=bins)
        pylab.legend()
        pylab.xlabel(xlabel, fontsize=fontsize)
        pylab.ylabel(ylabel, fontsize=fontsize)
        pylab.title(title,fontsize=fontsize)
        if grid is True:
            pylab.grid(True)

    def hist_nb_passes(self, bins=None, alpha=0.5, hold=False, fontsize=12,
                          grid=True, xlabel="Number of ZMW passes", logy=True,
                          ylabel="#", label="", title="Number of ZMW passes"):
        """Plot histogram of number of reads per ZMW (number of passes)

        :param float alpha: transparency of the histograms
        :param bool hold:
        :param int fontsize:
        :param bool grid:
        :param str xlabel:
        :param str ylabel:
        :param bool logy: use log scale on the y axis (default to True)
        :param str label: label of the histogram (for the legend)
        :param str title:

        .. plot::
            :include-source:

            from sequana.pacbio import PacbioSubreads
            from sequana import sequana_data
            b = PacbioSubreads(sequana_data("test_pacbio_subreads.bam"))
            b.hist_nb_passes()
        """
        max_nb_pass = self.df.nb_passes.max()
        if bins is None:
            k = range(1, max_nb_pass+1)

        # histogram nb passes
        if hold is False:
            pylab.clf()
        pylab.hist(self.df.nb_passes, bins=bins, alpha=alpha,
                   label=label, log=logy, width=1)
        if len(k) < 5:
            pylab.xticks(range(6), range(6))

        pylab.xlabel(xlabel, fontsize=fontsize)
        pylab.ylabel(ylabel, fontsize=fontsize)
        pylab.title(title, fontsize=fontsize)
        if grid is True:
            pylab.grid(True)

    def get_number_of_ccs(self, min_length=50, max_length=15000):
        query = "read_length>=@min_length and read_length <=@max_length"
        dd = self.df.query(query)
        return len(dd.ZMW.unique())

    def get_mean_nb_passes(self, min_length=50, max_length=15000):
        query = "read_length>=@min_length and read_length <=@max_length"
        dd = self.df.query(query).groupby("ZMW").agg(pylab.mean)["nb_passes"]
        return dd.mean()

    '''def hist_mean_ccs(self, bins=1000, **kwargs):
        """Group subreads by ZMW and plot mean of read length for each polymerase"""
        data = self.df[['read_length', 'ZMW']].groupby('ZMW')
        data.mean().hist(bins=bins, **kwargs)
        pylab.title("CCS mean read length")
        return data
    '''


class CCS(PacbioBAMBase):
    """

    You can get a CCS file from a BAM file as follows::

        singularity run /home/cokelaer/pacbio.simg ccs  --minPasses 80  \
            m54091_180224_083446.subreads.bam test.bam

    If empty, you may need to set other parameter such as the predicted
    accuracy::

    singularity run /home/cokelaer/pacbio.simg ccs --minPasses 0 \
        test_pacbio_subreads.bam out.bam  --minPredictedAccuracy 0.7

    """
    def __init__(self, filename):
        super(CCS, self).__init__(filename)

    @property
    def df(self):
        # RG: ID read group ??
        # np: number of passes
        # rq ?
        # rs: list 6 numbers ?
        # za:
        # zm ID of the ZMW
        # sn: SNR how is this computed ?
        # zs
        # - sn: list of ACGT SNRs. A, C, G, T in that order
        if self._df is not None:
            return self._df

        logger.info("Scanning input file. Please wait")
        self.reset()
        N = 0

        all_results = []
        # This takes 60%  of the time...could use cython ?
        for read in self.data:
            tags = dict(read.tags) #11% of the time
            res = []
            # count reads
            N += 1
            if (N % 10000) == 0:
                logger.info("Read %d sequences" %N)

            # res[0] = read length
            res.append(read.query_length) # also stored in tags["qe"] - tags["qs"]

            # collections.counter is slow, let us do it ourself
            res.append( 100. / read.qlen * sum(
                [read.query_sequence.count(letter) if read.query_sequence
                    else 0 for letter in "CGcgSs"]))

            # res[1:4] contains SNR  stored in tags['sn'] in the order A, C, G, T
            try:
                snr = list(tags['sn'])
            except:
                snr = [None] * 4
            res = res + snr

            # res[6] = ZMW name, also stored in tags["zm"]
            res.append(int(tags['zm']))
            res.append(tags['np'])

            # aggregate results
            all_results.append(res)

        self._df = pd.DataFrame(all_results,
            columns=['read_length','GC_content','snr_A','snr_C','snr_G','snr_T','ZMW',
                     "nb_passes"])
        self._df.ZMW = self._df.ZMW.astype(int)

        if len(self._df.ZMW.unique()) != len(self._df):
            logger.warning("Found non unique ZMW. This may not be a CCS but "
                        "a subread file. Consider using PacbioSubreads class")

        self.reset()
        return self._df

    def stats(self):
        data = {}
        data["N"] = len(self.df)
        data["mean_read_length"] = float(round(self.df.read_length.mean(),2))
        data["total_bases"] = int(self.df.read_length.sum())
        data["mean_nb_passes"] = float(round(self.df.nb_passes.mean(),2))
        return data

    def to_summary(self, filename="sequana_summary_pacbio_ccs.json"):
        """Save statistics into a JSON file

        :param filename:
        :param data: dictionary to save. If not provided, use :meth:`stats`

        """
        data = self.stats()
        s = Summary("pacbio_ccs", self.sample_name, data=data)
        s._data_description = {
            "N": "CCS reads",
            "total_bases": "Number of CCS bases",
            "mean_read_length": "CCS Read Length (mean)",
            "mean_nb_passes": "Number of Passes (mean)"
        }
        s.to_json(filename)

    def hist_passes(self, maxp=50, fontsize=16):
        passes = self.df.nb_passes.copy()
        passes.clip(upper=maxp).hist(bins=maxp)
        pylab.xlim([0, maxp])
        pylab.ylabel("# count", fontsize=fontsize)
        pylab.xlabel("Passes (max {})".format(maxp), fontsize=fontsize)

    def boxplot_read_length_vs_passes(self, nmax=20, ax=None, whis=1.5, widths=0.6):
        dd = self.df.query("nb_passes<=@nmax")[["nb_passes", "read_length"]]
        axes = dd.boxplot(
            by="nb_passes",
            notch=False, widths=widths,
            meanline=True, showmeans=True, ax=ax, whis=whis,
            flierprops=dict(markerfacecolor="blue", alpha=0.1, markersize=8),
            boxprops=dict(linewidth=1.5,  color="y"),
            medianprops=dict(color="k",linewidth=2.5, linestyle="-"),
            meanprops=dict(color="red", linestyle="--", linewidth=1.5),
            capprops=dict(color="black", linewidth=2),
            patch_artist=False)


class PacbioMappedBAM(PacbioBAMBase):

    def __init__(self, filename, method):
        super(PacbioMappedBAM, self).__init__(filename)
        assert method in ["bwa", "blasr", "minimap2"]
        self.method = method

    def _get_data(self):
        # return list of lists
        # each list is made of 3 values: mapq, length, concordance
        from sequana import Cigar
        data = []
        self.reset()
        count = 0

        for align in self.data:
            mapq = align.mapq
            length = align.rlen
            if self.method in ["blasr", "minimap2"]:
                this = Cigar(align.cigarstring).stats()
                S, D, I, M = this[4] , this[2] , this[1], this[0]
                concordance = 1 - (D+I+S)/(D + I + M + S)
            else:
                this = align.get_cigar_stats()[0]
                error = this[-1]  # suppose to be I + D + X
                total = this[-1] + this[0]
                if total:concordance = 1- (error)/(total)
                else:concordance = 0
            data.append([mapq, length, concordance])
            if count % 10000 == 0:
                logger.info("%s" % count)
            count+=1
        return data

    def filter_mapq(self, output_filename, threshold_min=0,
        threshold_max=255):
        """Select and Write reads within a given range

        :param str output_filename: name of output file
        :param int threshold_min: minimum length of the reads to keep
        :param int threshold_max: maximum length of the reads to keep

        """
        assert threshold_min < threshold_max
        assert output_filename != self.filename, \
            "output filename should be different from the input filename"
        self.reset()
        count = 0
        with pysam.AlignmentFile(output_filename,  "wb", template=self.data) as fh:
            for read in self.data:
                if ((read.mapq < threshold_max) & (read.mapq > threshold_min)):
                    fh.write(read)
                else:
                    pass
                count += 1
                if count % 10000:
                    logger.info("%s sequence processed" % count)

    def _set_concordance(self):
        from sequana import Cigar
        self._concordance = []
        self.reset()
        count = 0
        if self.method in ["blasr", "minimap2"]:
            for align in self.data:
                if align.cigarstring:
                    this = Cigar(align.cigarstring).as_dict()
                    S, D, I, M = this["S"] , this["D"] , this["I"], this["M"]
                    self._concordance.append((1- (D+I+S)/(D+I+M+S)))
                count += 1
                if count %1000 ==0 : print(count)
                if count == 10000: break
        elif self.method == "bwa":
            for align in self.data:
                if align.cigarstring:
                    this = align.get_cigar_stats()[0]
                    # Last item is total error though but cannot use Cigar class
                    # that does not have that value
                    error = this[-1]  # suppose to be I + D + X
                    # can check with alignment.get_cigar_stats() on blasr data
                    # for instance
                    #
                    total = this[-1] + this[0]
                    self._concordance.append((1- (error)/(total)))

    def hist_concordance(self,  bins=100, fontsize=16):
        """

            formula : 1 - (in + del + mismatch / (in + del + mismatch + match) )

        For BWA and BLASR, the get_cigar_stats are different !!!
        BWA for instance has no X stored while Pacbio forbids the use of the M
        (CMATCH) tag. Instead, it uses X (CDIFF) and = (CEQUAL) characters.

        Subread Accuracy: The post-mapping accuracy of the basecalls.
        Formula: [1 - (errors/subread length)], where errors = number of deletions +
        insertions + substitutions.

        """
        try:
            concordance = self._concordance
        except:
            self._set_concordance()
            concordance = self._concordance

        pylab.hist(concordance, bins=bins)
        pylab.grid()
        mu = np.mean(concordance)
        median = np.median(concordance)
        pylab.axvline(mu, color='r', alpha=0.5)
        pylab.axvline(median, color='r', alpha=0.5, ls="--")
        pylab.xlabel("concordance", fontsize=fontsize)

    def boxplot_mapq_concordance(self):
        # method can only be bwa for now
        assert self.method == "bwa"
        data = self._get_data()
        df = pd.DataFrame(data, columns=["mapq", "length", "concordance"])
        pylab.clf()
        pylab.boxplot([df[df.mapq == i]['concordance'] for i in range(1,61)])
        pylab.xlabel("mapq")
        pylab.ylabel("concordance")
        pylab.grid()
        tt = [10,20,30,40,50,60]
        pylab.xticks(tt, tt)

    def get_coverage(self, reference_length=None):
        self.reset()
        start = [this.reference_start for this in self.data]
        self.reset()
        end = [this.reference_end for this in self.data]
        if reference_length:
            N  = reference_length
        else:
            N = max([x for x in end if x])

        coverage = np.zeros(N)
        for x, y in zip(start, end):
            if y and x>=0 and y>=0: coverage[x:y] += 1
            else: pass
        return coverage

    def hist_median_ccs(self, bins=1000, **kwargs):
        """Group subreads by ZMW and plot median of read length for each polymerase"""
        data = self.df[['read_length', 'ZMW']].groupby('ZMW')
        data.median().hist(bins=bins, **kwargs)
        pylab.title("CCS median read length")
        return data


class BAMSimul(PacbioBAMBase):
    """BAM reader for Pacbio simulated reads (PBsim)

    A summary of the data is stored in the attribute :attr:`df`. It contains
    information such as the length of the reads, the ACGT content, the GC content.

    """
    def __init__(self, filename):
        """.. rubric:: Constructor

        :param str filename: filename of the input pacbio BAM file. The content
            of the BAM file is not the ouput of a mapper. Instead, it is the
            output of a Pacbio (Sequel) sequencing (e.g., subreads).
        """
        super(BAMSimul, self).__init__(filename)

    def _get_df(self):
        if self._df is None:
            self.reset()
            N = 0

            all_results = []
            for read in self.data:
                res = []
                # count reads
                N += 1
                if (N % 10000) == 0:
                    logger.info("Read %d sequences" %N)
                #res[0] = read length
                res.append(read.query_length)
                # res[1] = GC content
                c = collections.Counter(read.query_sequence)
                res.append( 100 * (c['g'] + c['G'] + c['c'] + c['C']) /
                            float(sum(c.values())) )

                # aggregate results
                all_results.append(res)


            self._df = pd.DataFrame(all_results,
                columns=['read_length','GC_content'])
            self.reset()
        return self._df
    df = property(_get_df)

    def filter_length(self, output_filename, threshold_min=0,
        threshold_max=np.inf):
        """Select and Write reads within a given range

        :param str output_filename: name of output file
        :param int threshold_min: minimum length of the reads to keep
        :param int threshold_max: maximum length of the reads to keep

        """
        assert threshold_min < threshold_max
        assert output_filename != self.filename, \
            "output filename should be different from the input filename"
        self.reset()
        with pysam.AlignmentFile(output_filename,  "wb", template=self.data) as fh:
            for read in self.data:
                if ((read.query_length > threshold_min) & (read.query_length < threshold_max)):
                    fh.write(read)

    def filter_bool(self, output_filename, mask):
        """Select and Write reads using a mask

        :param str output_filename: name of output file
        :param list list_bool: True to write read to output, False to ignore it

        """
        assert output_filename != self.filename, \
            "output filename should be different from the input filename"
        assert len(mask) == len(self), \
            "list of bool must be the same size as BAM file"
        self.reset()
        with pysam.AlignmentFile(output_filename,  "wb", template=self.data) as fh:
            for read, keep in zip(self.data, mask):
                if keep:
                    fh.write(read)


class PBSim(object):
    """Filter an input BAM (simulated with pbsim) so as so keep
    reads that fit a target distribution.

    This uses a MH algorithm behind the scene.

    ::

        ss = pacbio.PBSim("test10X.bam")
        clf();
        ss.run(bins=100, step=50)


    For example, to simulate data set, use::

        pbsim --data-type CLR --accuracy-min 0.85 --depth 20  \
            --length-mean 8000 --length-sd 800 reference.fasta --model_qc model_qc_clr

    The file model_qc_clr can be retrieved from the github here below.

    See https://github.com/pfaucon/PBSIM-PacBio-Simulator for details.

    We get a fastq file where simulated read sequences are randomly sampled from
    the reference sequence ("reference.fasta") and differences (errors) of the
    sampled reads are introduced.

    The Fastq can be converted to

    """
    def __init__(self, input_bam, simul_bam):
        self.bam = PacbioSubreads(input_bam)
        self.Nreads = len(self.bam)
        self.bam_simul = BAMSimul(simul_bam)

    def target_distribution(self, xprime):

        """The target distribution

        Compute histogram. Get X, Y.  Given xprime, interpolate to get yprime
        use e.g. np.interp

        """
        return np.interp(xprime, self.X[1:self.bins+1], self.Y)

    def run(self, bins=50, xmin=0, xmax=30000, step=1000, burn=1000,alpha=1,output_filename=None):
        # compute histogram of the input reads once for all to be used
        # in the target_distribution method
        self.bins = bins
        self.Y, self.X = np.histogram(self.bam.df.read_length, bins=bins, density=True)

        lengths = self.bam_simul.df.read_length.values
        self.tokeep = []
        vec = []
        x = self.bam.df.read_length.mean()
        for i in range(self.bam_simul.df.shape[0]):
            can = lengths[i]
            aprob = min([alpha,self.target_distribution(can)/self.target_distribution(x)])
            #acceptance probability
            u = pylab.uniform(0,1)
            if u < aprob:
                x = can
                vec.append(x)
                self.tokeep.append(True)
            else:
                self.tokeep.append(False)

        #plotting the results:
        #theoretical curve
        x = pylab.arange(xmin, xmax, step)
        y = self.target_distribution(x)
        pylab.subplot(211)
        pylab.title('Metropolis-Hastings')
        pylab.plot(vec)
        pylab.subplot(212)

        pylab.hist(vec[burn:], bins=bins, density=1)
        pylab.plot(x,y,'r-')
        pylab.ylabel('Frequency')
        pylab.xlabel('x')
        pylab.legend(('PDF','Samples'))

        if output_filename is not None:
            self.bam_simul.filter_bool(output_filename, self.tokeep)



class Barcoding():
    """

    Read as input a file created by smrtlink that stores statistics about each
    barcode. This is a simple CSV file with one line per barcode<


    """
    def __init__(self, filename):
        self.df = pd.read_csv(filename)
        self.df_not_barcoded = self.df[self.df['Barcode Index'] == 'None']
        self.df_barcoded = self.df[self.df['Barcode Index'] != 'None']

    def plot_polymerase_per_barcode(self, fontsize=12, unbarcoded=True):
        """Number Of Polymerase Reads Per Barcode"""
        PR = self.df_barcoded["Polymerase Reads"].sum()
        data = self.df_barcoded['Polymerase Reads'].sort_values(ascending=False).values
        pylab.plot([int(x) for x in range(1, len(data)+1)], data, label="barcodes")
        pylab.axhline(data.mean(), color="r", label="average")

        try:
            if unbarcoded is True:
                unbar = self.df_not_barcoded['Polymerase Reads'].iloc[0]
                pylab.axhline(unbar, color="k", ls="--",label="not barcoded")
        except:
            pass

        pylab.xlabel("Barcode Rank Order", fontsize=fontsize)
        pylab.ylabel("Counts of Reads", fontsize=fontsize)
        pylab.title("Total Polymerase count: {}".format(PR))
        pylab.legend()
        pylab.ylim(ymin=0)
        try:pylab.tight_layout()
        except:pass

    def hist_polymerase_per_barcode(self, bins=10, fontsize=12):
        """histogram of number of polymerase per barcode

        Cumulative histogram gives total number of polymerase reads

        """
        PR = self.df_barcoded["Polymerase Reads"].sum()
        self.df_barcoded['Polymerase Reads'].hist(bins=bins, ec="k", rwidth=0.8)
        pylab.title("Total Polymerase count: {}".format(PR))
        pylab.xlabel("Number of Polymerase Reads", fontsize=fontsize)
        pylab.ylabel("Number of Barcoded Samples", fontsize=fontsize)
        try:pylab.tight_layout()
        except:pass

    def hist_mean_polymerase_read_length(self, bins=10, fontsize=12):
        self.df_barcoded['Mean Read Length'].hist(bins=bins, ec="k", rwidth=0.8)
        pylab.xlabel("Mean Polymerase Read Length", fontsize=fontsize)
        pylab.ylabel("Number of Barcoded Samples", fontsize=fontsize)
        try:pylab.tight_layout()
        except:pass

    def plot_subreads_histogram(self, bins=10, fontsize=12):
        self.df_barcoded['Subreads'].hist(bins=bins, ec="k", rwidth=0.8)
        pylab.xlabel("Number of subreads", fontsize=fontsize)
        pylab.ylabel("Number of Barcoded Samples", fontsize=fontsize)
        try:pylab.tight_layout()
        except:pass

    def hist_quality_per_barcode(self, bins=10, fontsize=12):
        self.df_barcoded['Mean Barcode Quality'].hist(bins=bins, ec="k", rwidth=0.8)
        pylab.xlabel("Mean Barcode Quality", fontsize=fontsize)
        pylab.ylabel("Number of Barcoded Samples", fontsize=fontsize)
        try:pylab.tight_layout()
        except:pass

    def __str__(self):
        PR = self.df_barcoded['Polymerase Reads']
        MRL = self.df_barcoded['Mean Read Length']
        PR_NOBC = self.df_not_barcoded['Polymerase Reads'].values[0]

        txt = "{} unique barcodes\n".format(len(self.df_barcoded))
        txt += "{} barcoded PR reads\n".format(PR.sum())
        txt += "{} unbarcoded PR Reads\n".format(PR_NOBC)
        txt += "{} total PR reads\n".format(int(PR.sum()+PR_NOBC))
        txt += "{} mean PR Reads per bar code\n".format(int(PR.mean()))
        txt += "{} max PR reads per barcode\n".format(PR.max())
        txt += "{} min PR reads per barcode\n".format(PR.min())

        PASSES = self.df_barcoded['Subreads'].sum() / self.df_barcoded['Polymerase Reads'].sum()
        txt += "{} min number of passes\n".format(int(PASSES))

        M = int((PR * MRL).sum()/PR.sum())
        txt += "{} Mean read length\n".format(M)

        txt += "{} Mean Longest Subread Length\n".format(
            int(self.df_barcoded['Longest Subread Length'].mean()))
        return txt

    def plot_and_save_all(self, dpi=100, directory="."):

        def savefile(filename):
            outname = directory + os.sep +filename
            pylab.savefig(outname, dpi=dpi)

        pylab.clf()
        self.hist_polymerase_per_barcode()
        savefile("barcoding_hist_polymerase_per_barcode.png")

        pylab.clf()
        self.hist_quality_per_barcode()
        savefile("barcoding_hist_quality_per_barcode.png")

        pylab.clf()
        self.hist_mean_polymerase_read_length()
        savefile("barcoding_hist_mean_polymerase_read_length.png")

        pylab.clf()
        self.plot_polymerase_per_barcode()
        savefile("barcoding_polymerase_per_barcode.png")

        pylab.clf()
        self.plot_subreads_histogram()
        savefile("barcoding_subreads_histogram.png")

        print(self)

