#
#  This file is part of Sequana software
#
#  Copyright (c) 2016-2021 - Sequana Development Team
#
#  Distributed under the terms of the 3-clause BSD license.
#  The full license is in the LICENSE file, distributed with this software.
#
#  website: https://github.com/sequana/sequana
#  documentation: http://sequana.readthedocs.io
#
##############################################################################
"""
Python script to filter a VCF file
"""
import sys

import vcfpy

from sequana.lazy import pylab

import colorlog

logger = colorlog.getLogger(__name__)


__all__ = ["VCFBase", "strand_ratio", "compute_frequency", "compute_strand_balance"]


class VCFBase:
    """Base class for VCF files

    Read an existing file as follows::

        from sequana.vcf_filter import VCFBase
        v = VCFBase("filename.vcf")

    You can get the number of variants::

        len(v)

    the version and source of the VCF creator (if provided)::

        v.version, v.source

    and you can easily iterate through the variants::

        for variant in vcf:
            print(variant)

    note that if you iterate again, you will get nothing. You will need to
    rewind the cursor::

        vcf.rewind()

    You also get lots of extra information inherited from the vcf.Reader

    """

    def __init__(self, filename, verbose=True, **kwargs):
        """.. rubric:: constructor

        :param str filename: a vcf file.
        :param kwargs: any arguments accepted by vcf.Reader

        """
        self.filename = filename
        self.rewind()

        if verbose:
            print("Found VCF version {}".format(self.version))

    def rewind(self):
        """Rewind the reader"""
        self._vcf_reader = vcfpy.Reader.from_path(self.filename)

    def __len__(self):
        self.rewind()
        i = 0
        for line in self:
            i += 1
        self.rewind()
        return i

    def __iter__(self):
        return self

    def __next__(self):
        return next(self._vcf_reader)

    def _get_version(self):
        fileformat = [x for x in self._vcf_reader.header.lines if x.key == "fileformat"]
        if len(fileformat) == 1:
            fileformat = fileformat[0].value
        else:
            fileformat = "unknown"

        if fileformat == "VCFv4.1":
            return "4.1"
        elif fileformat == "VCFv4.2":
            return "4.2"
        else:
            return fileformat

    version = property(_get_version)

    def _get_source(self):
        fileformat = [x for x in self._vcf_reader.header.lines if x.key == "source"]
        if len(fileformat) == 1:
            fileformat = fileformat[0].value
        else:
            fileformat = "unknown"

    source = property(_get_source)

    def hist_qual(self, fontsize=16, bins=100):
        """

        This uses the QUAL information to be found in the VCF and should
        work for all VCF with version 4.1 (at least)

        """
        # TODO: could be moved to VCFBase
        self.rewind()
        data = [x.QUAL for x in self._vcf_reader]
        pylab.hist(data, bins=bins)
        pylab.grid(True)
        pylab.xlabel("Variant quality", fontsize=fontsize)


def strand_ratio(number1, number2):
    """Compute ratio between two number. Return result between [0:0.5]."""
    try:
        division = float(number1) / (number1 + number2)
        if division > 0.5:
            division = 1 - division
    except ZeroDivisionError:
        return 0
    return division


def compute_frequency(record):
    """Compute frequency of alternate allele.
        alt_freq = Count Alternate Allele / Depth

    :param record: variant record
    """
    try:
        info = record.info
    except:
        info = record.INFO

    alt_freq = [float(count) / info["DP"] for count in info["AO"]]
    return alt_freq


def compute_strand_balance(record):
    """Compute strand balance of alternate allele include in [0,0.5].
    strand_bal = Alt Forward / (Alt Forward + Alt Reverse)

    :param record: variant record


    FYI: in freebayes, the allele balance (reported under AB), strand
    bias counts (SRF, SRR, SAF, SAR) and bias estimate (SAP)
    can be used as well for filtering. Here, we use the strand balance
    computed as SAF / (SAF + SAR)


    """
    try:
        info = record.info
    except:
        info = record.INFO

    strand_bal = [strand_ratio(info["SAF"][i], info["SAR"][i]) for i in range(len(info["SAF"]))]

    return strand_bal
