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
"""Data analysis tool


.. autosummary::

    RunningMedian


"""
from bisect import bisect_left, insort
import numpy as np

import colorlog
logger = colorlog.getLogger(__name__)




# blist seems to be unstable on older systems/platforms so we use list by
# default for now. Be aware that on recent systems blist exhibits a log(W)
# complexity that is better than list complexity. Note, however that there is an
# overhead so the list is faster for W<20,000, which is the case in most
# applications.

def running_median(data, width, container=list):
    rm = RunningMedian(data, width, container=list)
    return rm.run()


class RunningMedian:
    """Running median (fast)

    This is an efficient implementation of running median, faster than SciPy
    implementation v0.17 and a skip list method.

    The main idea comes from a recipe posted in this website:
    http://code.activestate.com/recipes/576930/#c3 that uses a simple list
    as proposed in https://gist.github.com/f0k/2f8402e4dfb6974bfcf1 and was
    adapted to our needs included object oriented implementation.

    .. note:: a circular running median is implemented in :class:`sequana.bedtools.GenomeCov`

    ::

        from sequana.running_median import RunningMedian
        rm = RunningMedian(data, 101)
        results = rm.run()


    .. warning:: the first W/2 and last W/2 positions should be ignored
        since they do not use W values. In this implementation, the last
        W/2 values are currently set to zero.

    This shows how the results agree with scipy

    .. plot::
        :include-source:

        from pylab import *
        import scipy.signal
        from sequana.running_median import RunningMedian

        clf()
        x = randn(100)
        plot(x, 'k')
        plot(RunningMedian(x,9).run(), 'r', lw=4)
        plot(scipy.signal.medfilt(x, 9), 'go')
        grid()


    .. plot::
        :include-source:

        from sequana.running_median import RunningMedian
        from pylab import *

        N = 1000
        X = linspace(0, N-1, N)

        # Create some interesting data with SNP and longer over
        # represented section.
        data = 20 + randn(N) + sin(X*2*pi/1000.*5)
        data[300:350] += 10
        data[500:505] += 100
        data[700] = 1000

        plot(X, data, "k", label="data")
        rm = RunningMedian(data, 101)
        plot(X, rm.run(), 'r', label="median W=201")

        from sequana.stats import moving_average as ma
        plot(X[100:-100], ma(data, 201), 'g', label="mean W=201")
        grid()
        legend()
        ylim([10, 50])

    Note that for visualisation, we set the ylimits to 50 but the data at
    position 500 goes up to 120 and there is an large outlier (1000) at
    position 700 .

    We see that the median is less sensible to the outliers, as expected. The
    median is highly interesting for large outliers on short duration (e.g. here
    the peak at position 500) but is also less biases by larger regions.


    .. note:: The beginning and end of the running median are irrelevant.
        There are actually equal to the data in our implementation.

    .. note:: using blist instead of list is not always faster. It depends on
        the width of the window being used. list and blist are equivaltn for W
        below 20,000 (list is slightly faster). However, for large W, blist
        has an O(log(n)) complexity while list has a O(n) complexity

    """
    def __init__(self, data, width, container=list):
        """.. rubric:: constructor

        :param data: your data vector
        :param width: running window length
        :param container: a container (defaults to list). Could be a B-tree
            blist from the blist package but is 30% slower than a pure list
            for W < 20,000

        scipy in O(n)
        list in sqrt(n)
        blist in O(log(n))

        """
        if (width % 2) != 1:
            print("Warning[sequana]:: window length should be odd. Added +1.")
            width += 1

        self.container = container
        self.W = width
        self.data = data

    def __call__(self):
        return self.run()

    def run(self):

        # initialise with first W values and sort the values
        lc = self.container(self.data[:self.W])
        lc.sort()  # 

        mididx = (self.W - 1) // 2
        result = np.empty_like(self.data)

        # We initialise the first element
        idx = mididx
        result[idx] = lc[mididx]

        # We start at position W removing first element in lc
        # and adding a new one. We do not use enumerate since we do not
        # start at zero.
        for new_elem in self.data[self.W:]:
            old_elem = self.data[idx-mididx]
            del lc[bisect_left(lc, old_elem)]
            insort(lc, new_elem)
            idx += 1
            result[idx] = lc[mididx]

        # We decided to keep the first W/2 and last W/2 values as in the
        # original data. Ideally, they should not be used for post processing
        result[0:mididx] = self.data[0:mididx]
        result[-mididx:] = self.data[-mididx:]

        return result



