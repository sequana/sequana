#
#  This file is part of Sequana software
#
#  Copyright (c) 2016-2022 - Sequana Development Team
#
#  Distributed under the terms of the 3-clause BSD license.
#  The full license is in the LICENSE file, distributed with this software.
#
#  website: https://github.com/sequana/sequana
#  documentation: http://sequana.readthedocs.io
#
##############################################################################
"""Utilities for the genome coverage"""
import copy
import gc
import os
import random
import sys
from collections import Counter

import colorlog
from easydev import Progress, TempFile
from tqdm import tqdm

from sequana.errors import BadFileFormat
from sequana.genbank import GenBank
from sequana.gff3 import GFF3
from sequana.lazy import numpy as np
from sequana.lazy import pandas as pd
from sequana.lazy import pylab, pysam
from sequana.stats import evenness
from sequana.summary import Summary

logger = colorlog.getLogger(__name__)


__all__ = ["SequanaCoverage", "ChromosomeCov", "DoubleThresholds"]


class DoubleThresholds(object):
    """Simple structure to handle the double threshold for negative and
    positive sides

    Used by the :class:`SequanaCoverage` and related classes.

    ::

        dt = DoubleThresholds(-3, 4, 0.5, 0.5)

    This means the low threshold is -3 while the high threshold is 4. The two
    following values must be between 0 and 1 and are used to define the value
    of the double threshold set to half the value of the main threshold.

    Internally, the main thresholds are stored in the low and high attributes.
    The secondary thresholds are derived from the main thresholds and the
    two ratios. The ratios are named ldtr and hdtr for low double threshold
    ratio and high double threshold ratio. The secondary thresholds are
    denoted low2 and high2 and are update automatically if low, high, ldtr or
    hdtr are changed.

    """

    def __init__(self, low=-3, high=3, ldtr=0.5, hdtr=0.5):
        assert ldtr >= 0.0 and ldtr <= 1.0, "ldrt parameter (low double threshold ratio) must be in [0,1]"
        assert hdtr >= 0.0 and hdtr <= 1.0, "hdrt parameter (high double threshold ratio) must be in [0,1]"
        assert low < 0, "low threshold must be negative"
        assert high > 0, "high threshold must be positive"

        self._ldtr = ldtr
        self._hdtr = hdtr
        self._high = high
        self._low = low

    def _get_ldtr(self):
        return self._ldtr

    def _set_ldtr(self, ldtr):
        self._ldtr = ldtr
        self._low2 = self._low * self._ldtr

    ldtr = property(_get_ldtr, _set_ldtr)

    def _get_hdtr(self):
        return self._hdtr

    def _set_hdtr(self, hdtr):
        self._hdtr = hdtr
        self._high2 = self._high * self._hdtr

    hdtr = property(_get_hdtr, _set_hdtr)

    def _get_low(self):
        return self._low

    def _set_low(self, value):
        assert value < 0.0
        self._low = value
        self._low2 = self._low * self._ldtr

    low = property(_get_low, _set_low)

    def _get_high(self):
        return self._high

    def _set_high(self, value):
        assert value > 0.0
        self._high = value
        self._high2 = self._high * self._ldtr

    high = property(_get_high, _set_high)

    def _get_low2(self):
        return self._low * self._ldtr

    low2 = property(_get_low2)

    def _get_high2(self):
        return self._high * self._hdtr

    high2 = property(_get_high2)

    def get_args(self):
        return "%.2f,%.2f,%.2f,%.2f" % (self.low, self.high, self.ldtr, self.hdtr)

    def copy(self):
        thresholds = DoubleThresholds(self.low, self.high, self.ldtr, self.hdtr)
        return thresholds

    def __str__(self):
        txt = f"Low threshold: {self.low}\n"
        txt += f"High threshold: {self.high}\n"
        txt += f"double-low threshold: {self.low2}\n"
        txt += f"double-high threshold: {self.high2}"
        return txt


class SequanaCoverage(object):
    """Create a list of dataframe to hold data from a BED file generated with
    samtools depth.

    This class can be used to plot the coverage resulting from a mapping, which
    is stored in BED format. The BED file may contain several chromosomes.
    There are handled independently and accessible as a list of
    :class:`ChromosomeCov` instances.

    Example:

    .. plot::
        :include-source:

        from sequana import SequanaCoverage, sequana_data

        filename = sequana_data('JB409847.bed')
        reference = sequana_data("JB409847.fasta")

        gencov = SequanaCoverage(filename)

        # you can change the thresholds:
        gencov.thresholds.low = -4
        gencov.thresholds.high = 4
        #gencov.compute_gc_content(reference)

        gencov = SequanaCoverage(filename)
        for chrom in gencov:
            chrom.running_median(n=3001, circular=True)
            chrom.compute_zscore()
            chrom.plot_coverage()
        gencov[0].plot_coverage()

    Results are stored in a list of :class:`ChromosomeCov` named
    :attr:`chr_list`. For Prokaryotes and small genomes, this API
    is convenient but takes lots of memory for larger genomes.

    Computational time information: scanning 24,000,000 rows

        - constructor (scanning 40,000,000 rows): 45s
        - select contig of 24,000,000 rows: 1min20
        - running median: 16s
        - compute zscore: 9s
        - c.get_rois() ():

    """

    def __init__(
        self,
        input_filename,
        annotation_file=None,
        low_threshold=-4,
        high_threshold=4,
        ldtr=0.5,
        hdtr=0.5,
        force=False,
        chunksize=5000000,
        quiet_progress=False,
        chromosome_list=[],
        reference_file=None,
        gc_window_size=101,
    ):
        """.. rubric:: constructor

        :param str input_filename: the input data with results of a bedtools
            genomecov run. This is just a 3-column file. The first column is a
            string (chromosome), second column is the base postion and third
            is the coverage.
        :param str annotation_file: annotation file of your reference (GFF3/Genbank).
        :param float low_threshold: threshold used to identify under-covered
            genomic region of interest (ROI). Must be negative
        :param float high_threshold: threshold used to identify over-covered
            genomic region of interest (ROI). Must be positive
        :param float ldtr: fraction of the low_threshold to be used to define
            the intermediate threshold in the double threshold method. Must be
            between 0 and 1.
        :param float rdtr: fraction of the low_threshold to be used to define
            the intermediate threshold in the double threshold method. Must be
            between 0 and 1.
        :param chunksize: size of segments to analyse. If a chromosome is larger
            than the chunk size, it is split into N chunks. The segments are
            analysed indepdently and ROIs and summary joined together. Note that
            GC, plotting functionalities will only plot the first chunk.
        :param force: some constraints are set in the code to prevent unwanted
            memory issues with specific data sets of parameters. Currently,
            by default, (i) you cannot set the threshold below 2.5 (considered as
            noise).
        :param chromosome_list: list of chromosomes to consider (names). This is useful for
            very large input data files (hundreds million of lines) where each
            chromosome can be analysed one by one. Used by the sequana_coverage
            standalone. The only advantage is to speed up the constructor creation
            and could also be used by the Snakemake implementation.
        :param reference_file: if provided, computes the GC content
        :param int gc_window_size: size of the GC sliding window. (default 101)
        """
        # Keep various information as attributes

        self._circular = None
        self._feature_dict = None
        self._gc_window_size = gc_window_size
        self._annotation_file = None
        self._reference_file = reference_file
        self._window_size = None
        self._input_filename = input_filename

        # place holder for basic stats on all chromosomes.
        self._stats = {}
        self._basic_stats = {}
        self._rois = {}
        self._html_list = set()

        # setters
        if annotation_file:
            self.annotation_file = annotation_file

        self.chunksize = chunksize
        self.force = force
        self.quiet_progress = quiet_progress

        if high_threshold < 2.5 and self.force is False:
            raise ValueError("high threshold must be >=2.5")
        if low_threshold > -2.5 and self.force is False:
            raise ValueError("low threshold must be <=-2.5")

        #
        self.thresholds = DoubleThresholds(low_threshold, high_threshold, ldtr, hdtr)

        #

        if not input_filename.endswith(".bed"):
            raise Exception(("Input file must be a BED file " "(chromosome/position/coverage columns"))

        # scan BED files to set chromosome names
        self._scan_bed()
        if chromosome_list:
            for chrom in chromosome_list:
                if chrom not in self.chrom_names:
                    raise ValueError(
                        f"You provided an incorrect chromosome name ({chrom}). Please chose one amongst {self.chrom_names}"
                    )
            self.chrom_names = chromosome_list

    def __getitem__(self, index):
        chrom = self.chrom_names[index]
        return ChromosomeCov(self, chrom, self.thresholds, self.chunksize)

    def __iter__(self):
        for chrom in self.chrom_names:
            yield ChromosomeCov(self, chrom, self.thresholds, self.chunksize)

    def __len__(self):
        return len(self.chrom_names)

    @property
    def input_filename(self):
        return self._input_filename

    @property
    def reference_file(self):
        return self._reference_file

    @property
    def circular(self):
        """Get the circularity of chromosome(s). It must be a boolean."""
        return self._circular

    @circular.setter
    def circular(self, circular):
        if isinstance(circular, bool):
            self._circular = circular
        else:
            logger.error("TypeError: Circular must be a boolean. True if your " "genome is circular and False if not.")
            sys.exit(1)

    @property
    def feature_dict(self):
        """Get the features dictionary of the genbank."""
        return self._feature_dict

    @property
    def gc_window_size(self):
        """Get or set the window size to compute the GC content."""
        return self._gc_window_size

    @gc_window_size.setter
    def gc_window_size(self, n):
        if n % 2 == 0:
            logger.warning("Window size should be an odd number. Incrementing by one.")
            self._gc_window_size = n + 1
        else:
            self._gc_window_size = n

    @property
    def annotation_file(self):
        """Get or set the genbank filename to annotate ROI detected with
        :meth:`ChromosomeCov.get_roi`. Changing the genbank filename will
        configure the :attr:`SequanaCoverage.feature_dict`.
        """
        return self._annotation_file

    @annotation_file.setter
    def annotation_file(self, annotation_file):
        if os.path.isfile(annotation_file):
            self._annotation_file = os.path.realpath(annotation_file)

            try:
                gff = GFF3(annotation_file)
                self._feature_dict = gff.get_features_dict()
            except BadFileFormat:
                self._feature_dict = GenBank(annotation_file).genbank_features_parser()

                # just to be equivalent to the GFF, we fill genes with locus tag
                # See get_simplify_dataframe from GFF3 class. FIXME should be uniformed
                for (
                    chrom_name,
                    values,
                ) in self._feature_dict.items():  # loop over chromosomes
                    # values is a list of dictionaries
                    for i, value in enumerate(values):
                        if "gene" not in value and "locus_tag" in value:
                            self._feature_dict[chrom_name][i]["gene"] = value["locus_tag"]

        else:
            logger.error("FileNotFoundError: The genbank file doesn't exist.")
            sys.exit(1)

    @property
    def window_size(self):
        """Get or set the window size to compute the running median. Size
        must be an interger.
        """
        return self._window_size

    @window_size.setter
    def window_size(self, n):
        if n % 2 == 0:
            self._window_size = n + 1
        else:
            self._window_size = n

    def _scan_bed(self):
        # Figure out the length of the data and unique chromosome
        # name. This is required for large data sets. We do not store
        # any data here but the chromosome names and their positions
        # in the BED file. We also store the total length
        # Attributes filled:
        #  - chrom_names: list of contig/chromosome names
        #  - positions: dictionary with contig names and the starting and ending
        #               positions. Starting is zero in general, but not
        #               compulsary
        #  - total_length: number of rows in the BED file

        N = 0
        logger.info("Scanning input file (chunk of {} rows)".format(self.chunksize))
        positions = {}
        self.chrom_names = []

        # rough estimate of the number of rows for the progress bar
        logger.info("Estimating number of chunks to read")
        fullsize = os.path.getsize(self.input_filename)
        data = pd.read_table(self.input_filename, header=None, sep="\t", nrows=self.chunksize)
        with TempFile() as fh:
            data.to_csv(fh.name, index=None, sep="\t")
            smallsize = os.path.getsize(fh.name)
        del data

        Nchunk = int(fullsize / smallsize)
        i = 0
        for chunk in tqdm(
            pd.read_table(
                self.input_filename,
                header=None,
                sep="\t",
                usecols=[0],
                chunksize=self.chunksize,
                dtype="string",
            ),
            total=Nchunk,
            disable=self.quiet_progress,
        ):
            # accumulate length
            N += len(chunk)

            # keep unique names (ordered)
            contigs = chunk[0].unique()
            for contig in contigs:
                if contig not in self.chrom_names:
                    self.chrom_names.append(contig)

            # group by names (unordered)
            chunk = chunk.groupby(0)
            for contig in contigs:
                if contig not in positions:
                    positions[contig] = {
                        "start": chunk.groups[contig].min(),
                        "end": chunk.groups[contig].max(),
                    }
                else:
                    positions[contig]["end"] = chunk.groups[contig].max()
            i += 1
            i = min(i, Nchunk)
        del chunk
        self.total_length = N

        # Used in the main standalone
        for k in positions.keys():
            positions[k]["N"] = positions[k]["end"] - positions[k]["start"] + 1
        self.positions = positions

    def get_stats(self):
        """Return basic statistics for each chromosome

        :return: dictionary with chromosome names as keys
            and statistics as values.

        .. seealso:: :class:`ChromosomeCov`.

        .. note:: used in sequana_summary standalone
        """
        stats = {}
        for chrom_name in tqdm(self.chrom_names):
            if chrom_name in self._stats:
                stats[chrom_name] = self._stats[chrom_name]
            else:
                chrom = ChromosomeCov(self, chrom_name, self.thresholds, self.chunksize)
                stats[chrom.chrom_name] = chrom.get_stats()
                self._stats[chrom_name] = stats[chrom_name].copy()
        return stats


class ChromosomeCov(object):
    """Factory to manipulate coverage and extract region of interests.

    Example:

    .. plot::
        :include-source:

        from sequana import SequanaCoverage, sequana_data
        filename = sequana_data("virus.bed")

        gencov = SequanaCoverage(filename)

        chrcov = gencov[0]
        chrcov.running_median(n=3001)
        chrcov.compute_zscore()
        chrcov.plot_coverage()

        df = chrcov.get_rois().get_high_rois()

    The *df* variable contains a dataframe with high region of interests (over
    covered)


    If the data is large, the input data set is split into chunk. See
    :attr:`chunksize`, which is 5,000,000 by default.

    If your data is larger, then you should use the :meth:`run` method.

    .. seealso:: sequana_coverage standalone application
    """

    def __init__(self, genomecov, chrom_name, thresholds=None, chunksize=5000000):
        """.. rubric:: constructor

        :param df: dataframe with position for a chromosome used within
            :class:`SequanaCoverage`. Must contain the following columns:
            ["pos", "cov"]
        :param genomecov:
        :param chrom_name: to save space, no need to store the chrom name in the
            dataframe.
        :param thresholds: a data structure :class:`DoubleThresholds` that holds
            the double threshold values.
        :param chunksize: if the data is large, it is split and analysed by
            chunk. In such situations, you should use the :meth:`run` instead
            of calling the running_median and compute_zscore functions.

        """
        self._rois = None
        self.binning = 1

        # keep track of the chunksize user argument.
        self.chunksize = chunksize

        # and keep a link to the original SequanaCoverage instance.
        self._bed = genomecov

        # placeholder for GC vector
        self._gc_content = None

        # the chromosome of interest
        self.chrom_name = chrom_name

        # full data in memory or split into chunks ?
        self._mode = None

        self.thresholds = thresholds.copy()

        # open the file as an iterator (may be huge), set _df to None and
        # all attributes to None.
        self.init()

    def init(self):
        # jump to the file position corresponding to the chrom name.
        N = self.bed.positions[self.chrom_name]["N"]
        toskip = self.bed.positions[self.chrom_name]["start"]

        self.iterator = pd.read_table(
            self.bed.input_filename,
            skiprows=toskip,
            nrows=N,
            header=None,
            sep="\t",
            chunksize=self.chunksize,
            dtype={0: "string"},
        )

        if N <= self.chunksize:
            # we can load all data into memory:
            self._mode = "memory"
        else:
            self._mode = "chunks"

        self._df = None
        self._reset_metrics()

    def _reset_metrics(self):
        # attributes for stats
        self._evenness = None
        self._DOC = None
        self._BOC = None
        self._CV = None
        self._STD = None
        self._C3 = None
        self._C4 = None
        self._rois = None

    def __str__(self):
        N = self.bed.positions[self.chrom_name]["N"]
        if self._mode == "memory":
            stats = self.get_stats()
            txt = "\nChromosome length: {:>10}".format(N)
            txt += "\nSequencing depth (DOC):                {:>10.2f} ".format(self.DOC)
            txt += "\nSequencing depth (median):             {:>10.2f} ".format(stats["Median"])
            txt += "\nBreadth of coverage (BOC) (percent):   {:>10.2f} ".format(self.BOC)
            txt += "\nCoverage standard deviation:    {:>10.2f}".format(self.STD)
            txt += "\nCoverage coefficient variation: {:>10.2f}".format(self.CV)
        else:
            stats = self.get_stats()
            txt = "\nChromosome length: {:>10}".format(N)
            txt += "\n!!!! Information based on a sample of {} points".format(self.chunksize)
            txt += "\nSequencing depth (DOC):                {:>10.2f} ".format(self.DOC)
            txt += "\nSequencing depth (median):             {:>10.2f} ".format(stats["Median"])
            txt += "\nBreadth of coverage (BOC) (percent):   {:>10.2f} ".format(self.BOC)
            txt += "\nCoverage standard deviation:    {:>10.2f}".format(self.STD)
            txt += "\nCoverage coefficient variation: {:>10.2f}".format(self.CV)
            txt += "\nExact values will be available in the summary file (json)"
        return txt

    def __len__(self):

        # binning may be used, so len is not the length of the
        # dataframe, but the min/max positions.
        if self.binning > 1:
            return self._length
        else:
            return int(self.bed.positions[self.chrom_name]["N"])
        #    return self.df.__len__()

    def _set_chunk(self, chunk):
        chunk.rename(columns={0: "chr", 1: "pos", 2: "cov", 3: "mapq0"}, inplace=True)
        assert set(chunk["chr"].unique()) == set([self.chrom_name])
        chunk = chunk.set_index("chr", drop=True)
        chunk = chunk.set_index("pos", drop=False)
        self._df = chunk
        self._reset_metrics()
        # set GC if available
        if self.bed.reference_file:
            if self._gc_content is None:
                self._compute_gc_content()
            i1 = self._df.index[0] - 1
            i2 = self._df.index[-1]
            self._df["gc"] = self._gc_content[i1:i2]

    def next(self):
        try:
            chunk = next(self.iterator)
            self._set_chunk(chunk)
        except Exception as err:
            print(err)
            self._df = None

    def _check_window(self, W):
        if W * 2 > len(self.df):
            msg = "W ({}) is too large compared to the contig ({})".format(W, len(self.df))
            logger.error(msg)
            raise Exception(msg)

    def run(self, W, k=2, circular=False, binning=None, cnv_delta=None):

        self.init()
        # for the coverare snakemake pipeline
        if binning == -1:
            binning = None

        # for the coverare snakemake pipeline
        if cnv_delta == -1:
            cnv_delta = None

        if binning:
            logger.debug("binning={}".format(binning))
        if cnv_delta:
            logger.debug("cnv_delta={}".format(cnv_delta))
        if self.chunksize < 5 * W:
            logger.warning(("chunksize is small as compared to W. " "Not recommended! "))
        if binning:
            assert binning > 1 and isinstance(binning, int), "binning must be integer > 1"

        self.chunk_rois = []

        # Get the number of chunks
        num = self.bed.positions[self.chrom_name]["N"]
        denom = self.chunksize
        if num % denom == 0:
            N = num // denom
        else:
            N = num // denom + 1

        if binning is None or binning == 1:
            self.binning = 1
            if N > 1:
                pb = Progress(N)
                pb.animate(0)
            for i, chunk in enumerate(self.iterator):
                logger.debug("Analysing chunk {}".format(i + 1))
                self._set_chunk(chunk)

                logger.debug("running median computation")
                self.running_median(W, circular=circular)
                logger.debug("zscore computation")
                self.compute_zscore(k=k, verbose=False)  # avoid repetitive warning

                rois = self.get_rois()
                if cnv_delta is not None and cnv_delta > 1:
                    rois.merge_rois_into_cnvs(delta=cnv_delta)
                summary = self.get_summary()
                self.chunk_rois.append([summary, rois])
                if N > 1:
                    pb.animate(i + 1)
            if N > 1:
                print()
        else:
            # Need to store the length
            total_length = 0
            binned_df = None
            if N > 1:
                pb = Progress(N)
            for i, chunk in enumerate(self.iterator):
                # we group the data by small chunk of "binning" values
                self._set_chunk(chunk)
                NN = len(self._df)
                total_length += NN

                n = int(NN / binning)
                a = [i for i in range(n) for j in range(binning)]
                # if non-integer n, we need to fill the last ID manually
                # with the remaining points

                if len(a) < NN:
                    a = a + [n] * (NN - len(a))
                # Now, we groupby on the ID, and aggregate the cov columns (sum)
                # while the pos is aggregate by simply taking the minimum
                # (starting position).
                # self._set_chunk(chunk)
                df = self._df.copy()
                df["id"] = a
                groups = df.groupby("id")
                # Don't understand why but this is ultra slow
                # when called from this method, whereas, it is
                # fast in a shell, when calling statement by statement
                # Tried to use np.mean, pylab.mean, from np import mean,
                # the "mean" string command, etc...
                # df = groups.agg({
                #        "pos":lambda x: x.min(),
                #        "cov": np.mean
                #    })
                df = pd.DataFrame(
                    {
                        "pos": groups.min()["pos"].astype(int),
                        "cov": groups.mean()["cov"],
                    }
                )

                if "gc" in self._df.columns:
                    df["gc"] = groups.mean()["gc"]

                df.index = df["pos"]
                # append
                if binned_df is not None:
                    binned_df = pd.concat([binned_df, df])
                else:
                    binned_df = df
                if N > 1:
                    pb.animate(i + 1)
            if N > 1:
                print()
            # used by __len__
            self._length = total_length
            self._df = binned_df.copy()
            self.binning = binning

            self.running_median(int(W / binning), circular=circular)
            self.compute_zscore(k=k, verbose=False)  # avoid repetitive warning
            # Only one ROIs, but we use the same logic as in the chunk case,
            # and store the rois/summary in the ChromosomeCovMultiChunk
            # structure
            rois = self.get_rois()
            if cnv_delta is not None and cnv_delta > 1:
                rois.merge_rois_into_cnvs(delta=cnv_delta)
            summary = self.get_summary()
            self.chunk_rois.append([summary, rois])

        results = ChromosomeCovMultiChunk(self.chunk_rois)

        self._rois = results.get_rois()

        self.bed._basic_stats[self.chrom_name] = {
            "DOC": self.DOC,
            "CV": self.CV,
            "length": len(self),
        }
        # duplicated call to results.get_rois() FIXME
        self.bed._rois[self.chrom_name] = results.get_rois()

        # self._html_list = self._html_list.union(f"{sample}/{sample}.cov.html")

        return results

    @property
    def df(self):
        if self._df is None:
            self.next()
        return self._df

    @property
    def rois(self):
        if self._rois is None:
            self._rois = self.get_rois()
        return self._rois

    @property
    def bed(self):
        return self._bed

    def get_gaussians(self):
        return "{0}: {1}".format(self.chrom_name, self.gaussians_params)

    def moving_average(self, n, circular=False):
        """Compute moving average of the genome coverage

        :param n: window's size. Must be odd
        :param bool circular: is the chromosome circular or not

        Store the results in the :attr:`df` attribute (dataframe) with a
        column named *ma*.

        """
        N = len(self.df["cov"])
        assert n < N / 2
        from sequana.stats import moving_average

        ret = np.cumsum(np.array(self.df["cov"]), dtype=float)
        ret[n:] = ret[n:] - ret[:-n]
        ma = ret[n - 1 :] / n
        mid = int(n / 2)
        self._df["ma"] = pd.Series(ma, index=np.arange(start=mid, stop=(len(ma) + mid)))

        if circular:
            # FIXME: shift of +-1 as compared to non circular case...
            # shift the data and compute the moving average
            self.data = (
                list(self.df["cov"].values[N - n :]) + list(self.df["cov"].values) + list(self.df["cov"].values[0:n])
            )
            ma = moving_average(self.data, n)
            self.ma = ma[n // 2 + 1 : -n // 2]
            self._df["ma"] = pd.Series(self.ma, index=self.df["cov"].index)

    def running_median(self, n, circular=False):
        """Compute running median of genome coverage

        :param int n: window's size.
        :param bool circular: if a mapping is circular (e.g. bacteria
            whole genome sequencing), set to True

        Store the results in the :attr:`df` attribute (dataframe) with a
        column named *rm*.

        .. versionchanged:: 0.1.21
            Use Pandas rolling function to speed up computation.

        """
        self._check_window(n)
        self.circular = circular
        self.window_size = n

        # in py2/py3 the division (integer or not) has no impact
        mid = int(self.window_size / 2)
        self.range = [None, None]
        try:
            if self.circular:
                # BASED on running_median pure implementation, could be much
                # slower than pure pandas rolling function. Keep those 4 lines
                # for book keeping though.
                # cover = list(self.df["cov"])
                # cover = cover[-mid:] + cover + cover[:mid]
                # rm = running_median.RunningMedian(cover, n).run()
                # self.df["rm"] = rm[mid:-mid]
                rm = (
                    pd.concat([self.df["cov"][-mid:], self.df["cov"], self.df["cov"][:mid]])
                    .rolling(n, center=True)
                    .median()
                )
                self._df["rm"] = rm[mid:-mid]

            else:
                rm = self.df["cov"].rolling(n, center=True).median()
                # Like in RunningMedian, we copy the NAN with real data
                rm[0:mid] = self.df["cov"][0:mid]
                rm[-mid:] = self.df["cov"][-mid:]
                # rm = running_median.RunningMedian(cover, n).run()

                self._df["rm"] = rm
                # set up slice for gaussian prediction
                self.range = [mid, -mid]
        except:
            self._df["rm"] = self.df["cov"]

    @property
    def DOC(self):
        """depth of coverage"""
        if self._DOC is None:
            self._DOC = self.df["cov"].mean()
        return self._DOC

    @property
    def STD(self):
        """standard deviation of depth of coverage"""
        if self._STD is None:
            self._STD = self.df["cov"].std()
        return self._STD

    @property
    def CV(self):
        """The coefficient of variation (CV) is defined as sigma / mu"""
        if self._CV is None:
            if self.DOC == 0:
                self._CV = np.nan
            else:
                self._CV = self.STD / self.DOC
        return self._CV

    @property
    def BOC(self):
        """breadth of coverage"""
        if self._BOC is None:
            self._BOC = 100 * (1 - (self.df["cov"] == 0).sum() / float(len(self.df)))
        return self._BOC

    @property
    def C3(self):
        if self._C3 is None:
            self._C3 = self.get_centralness(3)
        return self._C3

    @property
    def C4(self):
        if self._C4 is None:
            self._C4 = self.get_centralness(4)
        return self._C4

    @property
    def evenness(self):
        """Return Evenness of the coverage

        :Reference: Konrad Oexle, Journal of Human Genetics 2016, Evaulation
            of the evenness score in NGS.

        work before or after normalisation but lead to different results.

        """
        if self._evenness is None:
            try:
                self._evenness = evenness(self.df["cov"])
            except:
                self._evenness = 0
        return self._evenness

    def _coverage_scaling(self):
        """Normalize data with moving average of coverage

        Store the results in the :attr:`df` attribute (dataframe) with a
        column named *scale*.

        .. note:: Needs to call :meth:`running_median`

        """
        if "rm" not in self.df.columns:
            txt = "Column rm (running median) is missing.\n" + self.__doc__
            print(txt)
            raise KeyError
        else:
            self._df["scale"] = self.df["cov"] / self.df["rm"]
        self._df = self.df.replace(np.inf, np.nan)
        self._df = self.df.replace(-np.inf, np.nan)

    def _get_best_gaussian(self):
        results_pis = [model["pi"] for model in self.gaussians_params]
        indice = np.argmax(results_pis)
        return self.gaussians_params[indice]

    def compute_zscore(self, k=2, use_em=True, clip=4, verbose=True):
        """Compute zscore of coverage and normalized coverage.

        :param int k: Number gaussian predicted in mixture (default = 2)
        :param float clip: ignore values above the clip threshold

        Store the results in the :attr:`df` attribute (dataframe) with a
        column named *zscore*.

        .. note:: needs to call :meth:`running_median` before hand.

        """
        # here for lazy import
        from sequana import mixture

        # normalize coverage
        self._coverage_scaling()

        # ignore start and end (corrupted due to running median window)
        # the slice does not seem to work as a copy, hence the copy()
        data = self.df["scale"][self.range[0] : self.range[1]].copy()

        # remove zero, nan and inf values and ignore values above 4 that would
        # bias the estimation of the central
        data[data > 4] = 0
        data = data.replace(0, np.nan)
        data = data.dropna()

        if data.empty:  # pragma: no cover
            self._df["scale"] = np.ones(len(self.df), dtype=int)
            self._df["zscore"] = np.zeros(len(self.df), dtype=int)
            # define arbitrary values
            self.gaussians_params = [
                {"mu": 0.5, "pi": 0.15, "sigma": 0.1},
                {"mu": 1, "pi": 0.85, "sigma": 0.1},
            ]
            self.best_gaussian = self._get_best_gaussian()

            return

        # if len data > 100,000 select 100,000 data points randomly
        if len(data) > 100000:
            indices = random.sample(range(len(data)), 100000)
            data = [data.iloc[i] for i in indices]

        if use_em:
            self.mixture_fitting = mixture.EM(data)
            self.mixture_fitting.estimate(k=k)
        else:
            self.mixture_fitting = mixture.GaussianMixtureFitting(data, k=k)
            self.mixture_fitting.estimate()

        # keep gaussians informations
        self.gaussians = self.mixture_fitting.results
        params_key = ("mus", "sigmas", "pis")
        self.gaussians_params = [{key[:-1]: self.gaussians[key][i] for key in params_key} for i in range(k)]
        self.best_gaussian = self._get_best_gaussian()

        # warning when sigma is equal to 0
        if self.best_gaussian["sigma"] == 0:
            logger.warning("A problem related to gaussian prediction is " "detected. Be careful, Sigma is equal to 0.")
            self._df["zscore"] = np.zeros(len(self.df), dtype=int)
        else:
            self._df["zscore"] = (self.df["scale"] - self.best_gaussian["mu"]) / self.best_gaussian["sigma"]

        # Naive checking that the 2 models are sensible (mus are different)
        if k == 2:
            mus = self.gaussians["mus"]
            sigmas = self.gaussians["sigmas"]

            index0 = mus.index(self.best_gaussian["mu"])
            if index0 == 0:
                mu1 = mus[1]
                s1 = sigmas[1]
                s0 = sigmas[0]
                mu0 = mus[0]
            else:
                mu1 = mus[0]
                s1 = sigmas[0]
                s0 = sigmas[1]
                mu0 = mus[1]
            if abs(mu0 - mu1) < s0:
                if verbose:
                    logger.warning(
                        (
                            "\nFound mixture model parameters (k=2) where "
                            " |mu0-mu1| < sigma0. k=1 could be a better choice."
                            "mu0={}, m1={}, sigma0={}, sigma1={}".format(mu0, mu1, s0, s1)
                        )
                    )

        # If coverage is zero, clearly, this is the lowest bound
        # Yet, for low depth of coverage (e.g. around 5), a deleted
        # gene may have a zscore between 0 and 3, which is considered noise
        # If we have deleted area, we want to make sure that the zscore is
        # large enough that it can be considered as an event.

        # Here, we set the zscore value to at least the value of the threshold

        self._df["zscore"] = self._df["zscore"].fillna(0)

        floor = self.df["cov"] == 0  # fastest than query, or ==
        newvalues = self.df.loc[floor, "zscore"].apply(lambda x: min(self.thresholds.low - 0.01, x))
        self._df.loc[floor, "zscore"] = newvalues

        # finally, since we compute the zscore, rois must be recomputed
        self._rois = None

    def get_centralness(self, threshold=3):
        """Proportion of central (normal) genome coverage

        This is 1 - (number of non normal data) / (total length)

        .. note:: depends on the thresholds attribute being used.
        .. note:: depends slightly on :math:`W` the running median window
        """
        # keep original thresholds
        temp_low = self.thresholds.low
        temp_high = self.thresholds.high

        self.thresholds.low = -threshold
        self.thresholds.high = threshold

        try:
            filtered = self.get_rois()
        except Exception as err:  # pragma: no cover
            raise (err)
        finally:
            # restore the original threshold
            self.thresholds.low = temp_low
            self.thresholds.high = temp_high
        Cplus = (filtered.get_high_rois()["size"]).sum()
        Cminus = (filtered.get_low_rois()["size"]).sum()
        Centralness = 1 - (Cplus + Cminus) / float(len(self))

        return Centralness

    def get_rois(self):
        """Keep positions with zscore outside of the thresholds range.

        :return: a dataframe from :class:`FilteredGenomeCov`

        .. note:: depends on the :attr:`thresholds` low and high values.
        """
        if "zscore" not in self.df.columns:
            logger.critical(
                (
                    "you must call running_median and compute_zscore."
                    "alternatively, the run() method does the two steps"
                    " at once"
                )
            )
            raise Exception
        features = self.bed.feature_dict

        try:
            second_high = self.thresholds.high2
            second_low = self.thresholds.low2
            query = "zscore > @second_high or zscore < @second_low"

            # in the genbank, the names appears as e.g. JB12345
            # but in the fasta or BED files, it may be something like
            # gi|269939526|emb|FN433596.1|
            # so they do not match. We can try to guess it
            alternative = None

            if features:
                if self.chrom_name not in features.keys():

                    alternative = [x for x in self.chrom_name.split("|") if x]
                    alternative = alternative[-1]  # assume the accession is last
                    alternative = alternative.split(".")[0]  # remove version
                    if alternative not in features.keys():
                        msg = """Chromosome name (%s) not found
                            in the genbank. Make sure the chromosome names in
                            the BAM/BED files are compatible with the genbank
                            content. Genbank files contains the following keys """
                        for this in features.keys():
                            msg += f"\n                        - {this}"
                        logger.warning(msg % self.chrom_name)

            data = self.df.query(query)
            data.insert(0, "chr", self.chrom_name)

            if features:
                if alternative:
                    return FilteredGenomeCov(data, self.thresholds, features[alternative], step=self.binning)
                else:
                    return FilteredGenomeCov(
                        data,
                        self.thresholds,
                        features[self.chrom_name],
                        step=self.binning,
                    )
            else:
                return FilteredGenomeCov(data, self.thresholds, step=self.binning)
        except KeyError:  # pragma: no cover
            logger.error(
                "Column zscore is missing in data frame.\n"
                "You must run compute_zscore before get low coverage."
                "\n\n",
                self.__doc__,
            )
            sys.exit(1)

    def plot_rois(
        self,
        x1,
        x2,
        set_ylimits=False,
        rois=None,
        fontsize=16,
        color_high="r",
        color_low="g",
        clf=True,
    ):
        if rois is None:
            rois = self.rois

        high = rois.get_high_rois().query("end>@x1 and start<@x2")
        low = rois.get_low_rois().query("end>@x1 and start<@x2")

        self.plot_coverage(
            x1=x1,
            x2=x2,
            set_ylimits=set_ylimits,
            sample=False,
            fontsize=fontsize,
            clf=clf,
        )

        for start, end, cov in zip(high.start, high.end, high.mean_cov):
            pylab.plot([start, end], [cov, cov], lw=2, color=color_high, marker="o")

        for start, end, cov in zip(low.start, low.end, low.mean_cov):
            pylab.plot([start, end], [cov, cov], lw=2, color=color_low, marker="o")
        pylab.xlim([x1, x2])

    def plot_coverage(
        self,
        filename=None,
        fontsize=16,
        rm_lw=1,
        rm_color="#0099cc",
        rm_label="Running median",
        th_lw=1,
        th_color="r",
        th_ls="--",
        main_color="k",
        main_lw=1,
        main_kwargs={},
        sample=True,
        set_ylimits=True,
        x1=None,
        x2=None,
        clf=True,
    ):
        """Plot coverage as a function of base position.

        :param filename:
        :param rm_lw: line width of the running median
        :param rm_color: line color of the running median
        :param rm_color: label for the running median
        :param th_lw: line width of the thresholds
        :param th_color: line color of the thresholds
        :param main_color: line color of the coverage
        :param main_lw: line width of the coverage
        :param sample: if there are more than 1 000 000 points, we
            use an integer step to skip data points. We can still plot
            all points at your own risk by setting this option to False

        :param set_ylimits: we want to focus on the "normal" coverage ignoring
            unsual excess. To do so, we set the yaxis range between 0 and a
            maximum value. This maximum value is set to the minimum between the
            10 times the mean coverage and 1.5 the maximum of the high coverage
            threshold curve. If you want to let the ylimits free, set this
            argument to False
        :param x1: restrict lower x value to x1
        :param x2: restrict lower x value to x2 (x2 must be greater than x1)

        .. note:: if there are more than 1,000,000 points, we show only
            1,000,000 by points. For instance for 5,000,000 points,

        In addition to the coverage, the running median and coverage confidence
        corresponding to the lower and upper  zscore thresholds are shown.

        .. note:: uses the thresholds attribute.
        """
        # z = (X/rm - \mu ) / sigma
        # some view to restrict number of points to look at:
        if x1 is None:
            x1 = self.df.iloc[0].name
        if x2 is None:
            x2 = self.df.iloc[-1].name
        assert x1 < x2

        df = self.df.loc[x1:x2]

        high_zcov = (self.thresholds.high * self.best_gaussian["sigma"] + self.best_gaussian["mu"]) * df["rm"]
        low_zcov = (self.thresholds.low * self.best_gaussian["sigma"] + self.best_gaussian["mu"]) * df["rm"]

        if clf:
            pylab.clf()
        ax = pylab.gca()
        ax.set_facecolor("#eeeeee")
        pylab.xlim(df["pos"].iloc[0], df["pos"].iloc[-1])
        axes = []
        labels = []

        # 1,000,000 points is already a lot for matplotlib. Let us restrict
        # the number of points to a million for now.
        if len(df) > 1000000 and sample is True:
            logger.info("sub sample data for plotting the coverage")
            NN = int(len(df) / 1000000)
        else:  # pragma: no cover
            NN = 1

        # the main coverage plot
        (p1,) = pylab.plot(
            df["cov"][::NN],
            color=main_color,
            label="Coverage",
            linewidth=main_lw,
            **main_kwargs,
        )
        axes.append(p1)
        labels.append("Coverage")

        # The running median plot
        if rm_lw > 0:
            (p2,) = pylab.plot(df["rm"][::NN], color=rm_color, linewidth=rm_lw, label=rm_label)
            axes.append(p2)
            labels.append(rm_label)

        # The threshold curves
        if th_lw > 0:
            (p3,) = pylab.plot(
                high_zcov[::NN],
                linewidth=th_lw,
                color=th_color,
                ls=th_ls,
                label="Thresholds",
            )
            (p4,) = pylab.plot(
                low_zcov[::NN],
                linewidth=th_lw,
                color=th_color,
                ls=th_ls,
                label="_nolegend_",
            )
            axes.append(p3)
            labels.append("Thresholds")

        pylab.legend(axes, labels, loc="best")
        pylab.xlabel("Position", fontsize=fontsize)
        pylab.ylabel("Per-base coverage", fontsize=fontsize)
        pylab.grid(True)

        # sometimes there are large coverage value that squeeze the plot.
        # Let us restrict it. We can say that 10 times the average coverage is
        # the maximum. We can save the original plot and the squeezed one.

        if set_ylimits is True:
            # m1 = high_zcov.max(skipna=True)
            m4 = high_zcov[high_zcov > 0]
            if len(m4) == 0:
                m4 = 3
            else:
                m4 = m4.max(skipna=True)
            # ignore values equal to zero to compute mean average
            m3 = df[df["cov"] > 0]["cov"].mean()

            pylab.ylim([0, min([m4 * 2, m3 * 10])])
        else:
            pylab.ylim([0, pylab.ylim()[1]])

        try:
            pylab.gcf().set_layout_engine("tight")
        except:  # pragma: no cover
            pass

        if filename:
            pylab.savefig(filename)

    def _set_bins(self, df, binwidth):
        try:
            bins = np.arange(min(df), max(df) + binwidth, binwidth)
        except ValueError:
            return 100
        if bins.any():
            return bins
        return 100

    def plot_hist_zscore(self, fontsize=16, filename=None, max_z=6, binwidth=0.5, **hist_kargs):
        """Barplot of the zscore values"""
        pylab.clf()
        bins = self._set_bins(self.df["zscore"][self.range[0] : self.range[1]], binwidth)
        self.df["zscore"][self.range[0] : self.range[1]].hist(grid=True, bins=bins, **hist_kargs)
        pylab.xlabel("Z-score", fontsize=fontsize)
        try:
            pylab.gcf().set_layout_engine("tight")
        except:  # pragma: no cover
            pass
        pylab.xlim([-max_z, max_z])
        if filename:
            pylab.savefig(filename)

    def plot_hist_normalized_coverage(self, filename=None, binwidth=0.05, max_z=3):
        """Barplot of the normalized coverage with gaussian fitting"""
        pylab.clf()
        # if there are a NaN -> can't set up binning
        d = self.df["scale"][self.range[0] : self.range[1]].dropna()
        # remove outlier -> plot crash if range between min and max is too high
        d = d[np.abs(d - d.mean()) <= (4 * d.std())]
        bins = self._set_bins(d, binwidth)
        try:
            self.mixture_fitting.data = d
            self.mixture_fitting.plot(self.gaussians_params, bins=bins, Xmin=0, Xmax=max_z)
        except (AttributeError, ZeroDivisionError):  # pragma: no cover
            return

        pylab.grid(True)
        pylab.xlim([0, max_z])
        pylab.xlabel("Normalised per-base coverage")

        try:
            pylab.gcf().set_layout_engine("tight")
        except:  # pragma: no cover
            pass
        if filename:
            pylab.savefig(filename)

    def plot_hist_coverage(
        self,
        logx=True,
        logy=True,
        fontsize=16,
        N=25,
        fignum=1,
        hold=False,
        alpha=0.8,
        ec="k",
        filename=None,
        zorder=10,
        **kw_hist,
    ):
        """

        :param N:
        :param ec:
        """
        if hold is False:
            pylab.figure(fignum)
            pylab.clf()
        ax = pylab.gca()
        ax.set_facecolor("#eeeeee")

        data = self.df["cov"].dropna().values

        maxcov = data.max()
        if logx is True and logy is True:
            bins = pylab.logspace(0, pylab.log10(maxcov), N)
            pylab.hist(
                data,
                bins=bins,
                log=True,
                label=self.chrom_name,
                alpha=alpha,
                ec=ec,
                zorder=zorder,
                **kw_hist,
            )
            pylab.semilogx()
            pylab.xlabel("Coverage (log scale)", fontsize=fontsize)
            pylab.ylabel("Count (log scale)", fontsize=fontsize)
        elif logx is False and logy is True:
            pylab.hist(
                data,
                bins=N,
                log=True,
                label=self.chrom_name,
                alpha=alpha,
                ec=ec,
                zorder=zorder,
                **kw_hist,
            )
            pylab.xlabel("Coverage", fontsize=fontsize)
            pylab.ylabel("Count (log scale)", fontsize=fontsize)
        elif logx is True and logy is False:
            bins = pylab.logspace(0, pylab.log10(maxcov), N)
            pylab.hist(
                data,
                bins=N,
                label=self.chrom_name,
                alpha=alpha,
                zorder=zorder,
                ec=ec,
                **kw_hist,
            )
            pylab.xlabel("Coverage (log scale)", fontsize=fontsize)
            pylab.ylabel("Count", fontsize=fontsize)
            pylab.semilogx()
        else:
            pylab.hist(
                data,
                bins=N,
                label=self.chrom_name,
                alpha=alpha,
                zorder=zorder,
                ec=ec,
                **kw_hist,
            )
            pylab.xlabel("Coverage", fontsize=fontsize)
            pylab.ylabel("Count", fontsize=fontsize)
        pylab.grid(True)
        if filename:
            pylab.savefig(filename)

    def to_csv(self, filename=None, start=None, stop=None, **kwargs):
        """Write CSV file of the dataframe.

        :param str filename: csv output filename. If None, return string.
        :param int start: start row index.
        :param int stop: stop row index.

        Params of :meth:`pandas.DataFrame.to_csv`:

        :param list columns: columns you want to write.
        :param bool header: determine if the header is written.
        :param bool index: determine if the index is written.
        :param str float_format: determine the float format.
        """
        # Create directory to avoid errno 2
        if filename:
            directory = os.path.dirname(os.path.realpath(filename))
            try:
                os.makedirs(directory)
            except FileExistsError:
                if os.path.isdir(directory):
                    pass
                else:
                    msg = f"{directory} exist and it is not a directory"
                    logger.error(msg)
                    raise FileExistsError
        return self.df.loc[start:stop].to_csv(filename, **kwargs)

    def plot_gc_vs_coverage(
        self,
        filename=None,
        bins=None,
        Nlevels=None,
        fontsize=20,
        norm="log",
        ymin=0,
        ymax=100,
        contour=True,
        cmap="BrBG",
        **kwargs,
    ):
        """Plot histogram 2D of the GC content versus coverage"""
        if Nlevels is None or Nlevels == 0:
            contour = False

        data = self.df[["cov", "gc"]].copy()
        data["gc"] *= 100
        data = data.dropna()
        if bins is None:
            bins = [
                100,
                min(
                    int(data["gc"].max() - data["gc"].min() + 1),
                    max(5, self.bed.gc_window_size - 10),
                ),
            ]
            bins[0] = max(10, min(bins[0], self.df["cov"].max()))

        # FIXME jan 2018 there is currently an annoying warning in Hist2D
        import warnings

        with warnings.catch_warnings():
            # warnings.simplefilter("ignore")
            from sequana.viz import Hist2D

            h2 = Hist2D(data)

            try:
                h2.plot(
                    bins=bins,
                    xlabel="Per-base coverage",
                    ylabel=r"GC content (%)",
                    Nlevels=Nlevels,
                    contour=contour,
                    norm=norm,
                    fontsize=fontsize,
                    cmap=cmap,
                    **kwargs,
                )
            except:  # pragma: no cover
                h2.plot(
                    bins=bins,
                    xlabel="Per-base coverage",
                    ylabel=r"GC content (%)",
                    cmap=cmap,
                    Nlevels=Nlevels,
                    contour=False,
                    norm=norm,
                    fontsize=fontsize,
                    **kwargs,
                )

        pylab.ylim([ymin, ymax])
        try:
            pylab.gcf().set_layout_engine("tight")
        except:  # pragma: no cover
            pass

        if filename:
            pylab.savefig(filename)

    def get_gc_correlation(self):
        """Return the correlation between the coverage and GC content

        The GC content is the one computed in :meth:`SequanaCoverage.compute_gc_content`
        (default window size is 101)

        """
        self.df["gc"] = self._gc_content[self.df.index[0] - 1 : self.df.index[-1]]
        C = self.df[["cov", "gc"]].corr().iloc[0, 1]
        del self.df["gc"]
        gc.collect()
        return C

    def _get_hist_data(self, bins=30):
        data = self.df["cov"].dropna()
        m = data.quantile(0.01)
        M = data.quantile(0.99)
        # we
        step = 1
        bins = pylab.arange(m, M, step)

        while len(bins) > 150:
            step *= 2
            bins = pylab.arange(m, M, step)

        try:
            (
                Y,
                X,
            ) = np.histogram(data, bins=bins, density=True)
            return {"X": list(X[1:]), "Y": list(Y)}
        except:  # pragma: no cover
            return {"X": [], "Y": []}

    def get_summary(self, C3=None, C4=None, stats=None, caller="sequana.bedtools"):
        if stats is None:
            stats = self.get_stats()

        ROI = len(self.rois)
        ROIlow = len(self.rois.get_low_rois())
        ROIhigh = len(self.rois.get_high_rois())

        fit = self._get_best_gaussian()

        d = {
            "evenness": round(self.evenness, 4),
            "C3": round(self.C3, 4),
            "C4": round(self.C4, 4),
            "BOC": self.BOC,
            "length": len(self),
            "DOC": self.DOC,
            "CV": self.CV,
            "chrom_name": self.chrom_name,
            "ROI": ROI,
            "ROI(low)": ROIlow,
            "ROI(high)": ROIhigh,
            "fit_mu": fit["mu"],
            "fit_sigma": fit["sigma"],
            "fit_pi": fit["pi"],
            "hist_coverage": self._get_hist_data(),
        }
        if "gc" in stats:
            d["GC"] = stats["gc"]

        # sample name will be the filename
        # chrom name is the chromosome or contig name
        # Fixes v0.8.0 get rid of the .bed extension
        sample_name = os.path.basename(self._bed.input_filename.replace(".bed", ""))
        summary = Summary("coverage", sample_name=sample_name, data=d, caller=caller)

        summary.data_description = {
            "BOC": "Breadth of Coverage",
            "DOC": "Depth of Coverage",
            "chrom_name": "name of the contig/chromosome",
            "CV": "coefficient of variation of the DOC",
            "ROI": "number of regions of interest found",
            "C3": "Centralness (1 - ratio outliers by genome length) using zscore of 3",
            "C4": "Centralness (1 - ratio outliers by genome length) using zscore of 4",
            "length": "length of the chromosome",
        }
        if "gc" in stats:
            summary.data_description["GC"] = "GC content in %"

        return summary

    def get_stats(self):
        """Return basic stats about the coverage data

        only "cov" column is required.

        :return: dictionary
        """
        MedianCOV = self.df["cov"].median()

        stats = {
            "DOC": self.DOC,  # depth of coverage
            "STD": self.STD,  # sigma of the DOC
            "Median": MedianCOV,  # median of the DOC
            "BOC": self.BOC,  # breadth of coverage
        }

        # coefficient of variation
        stats["CV"] = self.CV

        # median of the absolute median deviation
        stats["MAD"] = np.median(abs(MedianCOV - self.df["cov"]).dropna())

        # GC content
        if "gc" in self.df.columns:
            stats["GC"] = self.df["gc"].mean() * 100

        stats["length"] = len(self)

        return stats

    def _compute_gc_content(self, ws=None):
        if self.bed.reference_file is None:
            return

        fasta = pysam.FastxFile(self.bed.reference_file)
        chrom_gc_content = dict()

        if ws is None:
            gc_window_size = self.bed.gc_window_size
        else:
            gc_window_size = ws

        mid = int(gc_window_size / 2.0)

        for chrom in fasta:
            if chrom.name == self.chrom_name:
                # Create gc_content array
                gc_content = np.empty(len(chrom.sequence))
                gc_content[:] = np.nan
                if self.bed.circular:
                    chrom.sequence = chrom.sequence[-mid:] + chrom.sequence + chrom.sequence[:mid]
                    # Does not shift index of array
                    mid = 0

                # Count first window content
                counter = Counter(chrom.sequence[0:gc_window_size])
                gc_count = 0
                for letter in "GCgc":
                    gc_count += counter[letter]

                gc_content[mid] = gc_count

                for i in range(1, len(chrom.sequence) - gc_window_size + 1):
                    if chrom.sequence[i - 1] in "GCgc":
                        gc_count -= 1
                    if chrom.sequence[i + gc_window_size - 1] in "GCgc":
                        gc_count += 1
                    gc_content[i + mid] = gc_count
                chrom_gc_content[chrom.name] = gc_content / gc_window_size

        self._gc_content = chrom_gc_content[self.chrom_name]


class FilteredGenomeCov(object):
    """Select a subset of :class:`ChromosomeCov`

    :target: developers only
    """

    _feature_not_wanted = {"regulatory", "source"}

    def __init__(
        self,
        df,
        threshold,
        feature_list=None,
        step=1,
        apply_threshold_after_merging=True,
    ):
        """.. rubric:: constructor

        :param df: dataframe with filtered position used within
            :class:`SequanaCoverage`. Must contain the following columns:
            ["pos", "cov", "rm", "zscore"]
        :param int threshold: a :class:`~sequana.bedtools.DoubleThresholds`
            instance.
        :param apply_threshold_after_merging: see :meth:`merge_rois_into_cnvs`.


        """
        self.rawdf = df.copy()
        self.rawdf
        self.thresholds = threshold
        self.apply_threshold_after_merging = True

        if isinstance(feature_list, list) and len(feature_list) == 0:
            feature_list = None

        self.feature_list = feature_list

        self.step = step
        region_list = self._merge_region()

        if self.feature_list:
            region_list = self._add_annotation(region_list, self.feature_list)

        self.df = self._dict_to_df(region_list, self.feature_list)

        def func(x):
            try:
                return x.split(".")[0]
            except:
                return x

        for column in ["gene_end", "gene_start"]:
            if column in self.df.columns:
                self.df[column] = self.df[column].astype(str)
                self.df[column] = self.df[column].apply(func)

    def __iadd__(self, other):
        self.df = self.df.append(other.df)
        return self

    def __str__(self):
        return self.df.__str__()

    def __len__(self):
        return self.df.__len__()

    def _merge_row(self, start, stop, chrom=None):
        df = self.rawdf

        chrom = df["chr"][start]

        cov = np.mean(df["cov"].loc[start:stop])
        max_cov = np.max(df["cov"].loc[start:stop])
        rm = np.mean(df["rm"].loc[start:stop])
        zscore = np.mean(df["zscore"].loc[start:stop])
        if zscore >= 0:
            max_zscore = df["zscore"].loc[start:stop].max()
        else:
            max_zscore = df["zscore"].loc[start:stop].min()
        size = stop - start + 1

        def log2ratio(mean_cov, rm):
            if rm != 0 and mean_cov != 0:
                return np.log2(mean_cov / rm)
            else:
                np.nan

        return {
            "chr": chrom,
            "start": start,
            "end": stop + 1,
            "size": size,
            "mean_cov": cov,
            "mean_rm": rm,
            "mean_zscore": zscore,
            "log2_ratio": log2ratio(cov, rm),
            "max_zscore": max_zscore,
            "max_cov": max_cov,
        }

    def _merge_region(self, zscore_label="zscore"):
        """Cluster regions within a dataframe.

        Uses a double thresholds method using the :attr:`threshold`

        """
        region_start = None
        region_stop = None
        start = 1
        stop = 1
        prev = 1
        # handle case where for example position n-1 have a zscore of -5 and n
        # have a zscore of 5. It is two different regions.
        region_zscore = 0

        merge_df = []
        for pos, zscore in zip(self.rawdf["pos"], self.rawdf[zscore_label]):
            stop = pos
            if stop - self.step == prev and zscore * region_zscore >= 0:
                prev = stop
            else:
                if region_start:
                    merge_df.append(self._merge_row(region_start, region_stop))
                    region_start = None
                start = stop
                prev = stop
                region_zscore = zscore

            if zscore > 0 and zscore > self.thresholds.high:
                if not region_start:
                    region_start = pos
                    region_stop = pos
                else:
                    region_stop = pos
            elif zscore < 0 and zscore < self.thresholds.low:
                if not region_start:
                    region_start = pos
                    region_stop = pos
                else:
                    region_stop = pos

        if start < stop and region_start:
            merge_df.append(self._merge_row(region_start, region_stop))
        return merge_df

    def _add_annotation(self, region_list, feature_list):
        """Add annotation from a dictionary generated by parsers in
        sequana.tools.
        """
        region_ann = []
        # an iterator of features
        iter_feature = iter(feature_list)
        feature = next(iter_feature)
        # pass "source" feature
        while feature["type"] in FilteredGenomeCov._feature_not_wanted:
            try:
                feature = next(iter_feature)
            except StopIteration:  # pragma: no cover
                print(
                    "Features types ({0}) are not present in the annotation"
                    " file. Please change what types you want".format(feature["type"])
                )
                return region_ann

        # merge regions and annotations
        for region in region_list:
            feature_exist = False
            while feature["gene_end"] <= region["start"]:
                try:
                    feature = next(iter_feature)
                except StopIteration:
                    break
            while feature["gene_start"] < region["end"]:
                # A feature exist for detected ROI
                feature_exist = True
                for item in ["gene", "product", "gene_name", "gene_id"]:
                    if item not in feature:
                        feature[item] = "not avail."

                region_ann.append(dict(region, **feature))
                try:
                    feature = next(iter_feature)
                except StopIteration:
                    break

            if not feature_exist:
                region_ann.append(
                    dict(
                        region,
                        **{
                            "gene_start": None,
                            "gene_end": None,
                            "type": None,
                            "gene": None,
                            "strand": None,
                            "product": None,
                        },
                    )
                )
        return region_ann

    def _dict_to_df(self, region_list, annotation):
        """Convert dictionary as dataframe."""
        merge_df = pd.DataFrame(region_list)
        colnames = [
            "chr",
            "start",
            "end",
            "size",
            "mean_cov",
            "max_cov",
            "mean_rm",
            "mean_zscore",
            "max_zscore",
            "log2_ratio",
            "gene_start",
            "gene_end",
            "type",
            "gene",
            "gene_name",
            "gene_id",
            "strand",
            "product",
        ]
        if not annotation:
            colnames = colnames[:10]
        merge_df = pd.DataFrame(region_list, columns=colnames)
        int_column = ["start", "end", "size"]
        merge_df[int_column] = merge_df[int_column].astype(int)
        if annotation:
            # maybe let the user set what he wants
            return merge_df.loc[~merge_df["type"].isin(FilteredGenomeCov._feature_not_wanted)]
        return merge_df

    def _get_sub_range(self, seq_range):
        try:
            return self.df[(self.df["end"] > seq_range[0]) & (self.df["start"] < seq_range[1])]
        except TypeError:
            return self.df

    def get_low_rois(self, seq_range=None):
        df = self._get_sub_range(seq_range)
        return df.loc[df["max_zscore"] < 0]

    def get_high_rois(self, seq_range=None):
        df = self._get_sub_range(seq_range)
        return df.loc[df["max_zscore"] >= 0]

    def merge_rois_into_cnvs(self, delta=1000):
        """

        E-E-E---------E-------E-E

        should be clustered as:
        <EEE>---------E-------<EE>

        Where new large events will combine information from
        the small events.

        When merging two events, we recompute the  max_zscore and mean_zscore
        and apply the threshold again.

        To illustrate this feature, imagine two marginal events at position 1
        and 1000 with a zscore of 4 each (marginals). When there are merged,
        their mean_zscore is recomputed and should be smaller than 4. Those
         events are just noise but will appear as significant events (since
        long size). So, we can remove then by checking the new mean_zscore to
        be > than the thresholds.

        """
        logger.info("Merging ROIs into CNV-like events")

        # merge using "CNV" clustering
        rois = self.df.copy()
        logger.info("Found {} ROIs".format(len(rois)))
        data = self._merge_rois_into_cnvs(rois, delta=delta)
        logger.info("Merged into {} ROIs".format(len(data)))

        # sort data by chr and start
        self.df = data.sort_values(by=["chr", "start"])

        # some rows have been merged, let us reset the index
        self.df.reset_index(inplace=True, drop=True)

    def _merge_rois_into_cnvs(self, rois, delta=1000):
        # rois is a copy, it can be changed

        clusterID = [-1] * len(rois)
        cluster_counter = 0
        for i in range(0, len(rois) - 2):
            start = rois.iloc[i : i + 2]["start"].values
            end = rois.iloc[i : i + 2]["end"].values

            # for the next ROI to be added as an item of
            # the cluster, it should be close, and similar.
            # similar for now means zscore have the same sign
            d1 = start[1] - end[0]
            z1 = rois.iloc[i]["max_zscore"]
            z2 = rois.iloc[i + 1]["max_zscore"]

            if d1 < delta and z1 * z2 >= 0:  # and d2 <1000:
                clusterID[i] = cluster_counter
                clusterID[i + 1] = cluster_counter
            else:
                cluster_counter += 1

        # Now, we can cluster the events
        # newdata contains the unclustered rows, then
        # we will add one row per cluster.
        rois["cluster"] = clusterID
        self.clusterID = clusterID
        self.rois = rois

        # just get start and end time of the clusters.
        # we drop special case of unclustered rows (-1)
        new_regions = rois.groupby("cluster").agg(
            {
                # min of ("A", "A") should return "A"
                # "chr": lambda x: x.min(),
                "start": lambda x: x.min(),
                "end": lambda x: x.max(),
            }
        )
        new_regions.drop(-1, inplace=True)

        # Now we build back the entire ROI dataframe
        region_list = []

        for _, row in new_regions.iterrows():
            this_cluster = self._merge_row(row.start, row.end)
            region_list.append(this_cluster)

        for _, row in rois.query("cluster == -1").iterrows():
            this_cluster = self._merge_row(row.start, row.end)
            region_list.append(this_cluster)

        merge_df = self._dict_to_df(region_list, self.feature_list)

        # finally, remove events that are small.
        if self.apply_threshold_after_merging:
            merge_df = merge_df.query("mean_zscore>@self.thresholds.high or mean_zscore<@self.thresholds.low")

        return merge_df


class ChromosomeCovMultiChunk(object):
    """Contains several chunk of :class:`ChromosomeCov` results.


    If genome/chromosome are too large, they cannot fit into memory.
    We therefore implemented a way to analyse the chromosome in chunks.

    This is done in :class:`ChromosomeCov`. Results, which are summaries and
    ROIs, are then stored into lists. For each chunk correspond a summary and a
    list of Regions of Interests. This is not convenient and needs extra
    processing.

    For instance, to get the final depth of coverage, (DOC), we need to retrieve
    the DOC for each chunk.

    Individual results are stored in :attr:`data` as a list of lists where each
    item contains one summary, and one data structure for the ROI.

    To retrieve all summaries, one would iterate as follows through the data::

        cc = ChromosomeCovMultiChunk(data)
        summaries = [item[0] for item in cc]

    and to retrieve all ROIs::

        summaries = [item[1] for item in cc]

    but one really wants to extra from this data structure is the overall
    summary::

        cc.get_summary()

    and concatenated list of ROIs::

        cc.get_rois()


    """

    def __init__(self, chunk_rois):
        self.data = chunk_rois

    def get_summary(self, caller="sequana.bedtools"):
        # get all summaries
        summaries = [this[0].as_dict() for this in self.data]

        # from the first one extract metadata
        data = summaries[0]
        sample_name = data["sample_name"]
        summary = Summary("coverage", sample_name=sample_name, data=data["data"].copy(), caller=caller)
        summary.data_description = data["data_description"].copy()

        # now replace the data field with proper values
        N = sum([this["data"]["length"] for this in summaries])
        summary.data["length"] = N

        # now, we need to update those values, which are means, so
        # we need to multiply back by the length to get the sum, and finally
        # divide by the total mean
        for this in ["BOC", "DOC", "evenness"]:
            summary.data[this] = sum([d["data"][this] * d["data"]["length"] for d in summaries]) / float(N)

        # For, CV, centralness, evenness, we simply takes the grand mean for now
        for this in ["C3", "C4", "evenness"]:
            summary.data[this] = np.mean([d["data"][this] for d in summaries])

        # For, ROI, just the sum
        for this in ["ROI", "ROI(high)", "ROI(low)"]:
            summary.data[this] = sum([d["data"][this] for d in summaries])

        return summary

    def get_rois(self):

        # all individual ROIs
        data = [item[1] for item in self.data]

        # let us copy the first one
        rois = copy.deepcopy(data[0])
        rois.df = pd.concat([this.df for this in data], ignore_index=True)
        return rois
