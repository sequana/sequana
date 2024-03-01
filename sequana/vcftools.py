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

import colorlog

from sequana.lazy import pylab
from sequana.utils.fisher import fisher_exact

logger = colorlog.getLogger(__name__)


__all__ = [
    "strand_ratio",
    "compute_frequency",
    "compute_strand_balance",
    "compute_fisher_strand_filter",
]


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

    return [float(count) / info["DP"] for count in info["AO"]]


def compute_fisher_strand_filter(record):
    """Fisher strand filter (FS).

    Sites where the numbers of reference/non-reference reads are highly correlated
    with the strands of the reads. Counting the number of reference reads on the forward
    strand and on the reverse strand, and the number of alternate reads on the forward and
    reverse strand should be equivalent. With these four numbers, we
    construct a 2 x 2 contingency table and used the P-value from a Fisher’s exact test
    to evaluate the correlation.

        from sequana import freebayes_vcf_filter
        v = freebayes_vcf_filter.VCF_freebayes("data/JB409847.vcf")
        r = next(v)
        compute_fisher_strand_filter(r)

    If the pvalue is less than 0.05 we should reject the variant since the alternate and reference
    do not behave in the same way. Typically found if the alternate has a poor strand balance. Used with
    care if frequency of alternate or reference is low or deth of coverage is low.


    .. note:: fisher terst with two-sided hypothesis"""
    try:
        info = record.info
    except:
        info = record.INFO

    return [fisher_exact([[x, y], [info["SRF"], info["SRR"]]], "two-sided") for x, y in zip(info["SAF"], info["SAR"])]


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

    return [strand_ratio(x, y) for x, y in zip(info["SAF"], info["SAR"])]
