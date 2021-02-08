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
"""Statistical tools"""

from sequana.lazy import numpy as np
from sequana.lazy import pandas as pd

import colorlog
logger = colorlog.getLogger(__name__)


__all__ = ["moving_average", "evenness", "runmean"]


def runmean(data, n):
    return moving_average(data, n)

def moving_average(data, n):
    """Compute moving average

    :param n: window's size (odd or even).

    ::

        >>> from sequana.stats import moving_average as ma
        >>> ma([1,1,1,1,3,3,3,3], 4)
        array([ 1. ,  1.5,  2. ,  2.5,  3. ])

    .. note:: the final vector does not have the same size as the input
        vector.

    """
    ret = np.cumsum(data, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    ma = ret[n - 1:] / n
    return ma


def evenness(data):
    """Return Evenness of the coverage

    :Reference: Konrad Oexle, Journal of Human Genetics 2016, Evaulation
        of the evenness score in NGS.

    work before or after normalisation but lead to different results.

    .. math::

        C = mean(X)
        D2 = X[X<=C]
        N = len(X)
        n = len(D2)
        E = 1 - (n - sum(D2) / C) / N

    """
    coverage = pd.Series(data)

    coverage = coverage.dropna()

    C = float(round(coverage.mean()))
    D2 = coverage[coverage<=C]
    if len(D2) == 0: #pragma: no cover
        return 1
    else:

        return 1. - (len(D2) - sum(D2) / C) / len(coverage)


def N50(data):
    """Return the N50 value given a list of unsorted/sorted contigs

    Once the list of contigs is sorted, the N50 is the contig length for which at
    least half of the nucleotides in the assembly belongs to contigs with the N50
    length or longer.

    """
    data = np.sort(data)
    cdata = np.cumsum(data)
    return data[np.argmax(cdata>cdata[-1]/2)]


def L50(data):
    """Return the smallest number of contigs whose length sum produces N50

    ::

        >>> data =
        >>> L50(data)
        3
    """
    data = np.sort(data)
    cdata = np.cumsum(data)
    pos = np.argmax(cdata>cdata[-1]/2)
    return len(data) - pos
