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
"""2D histograms



"""
from sequana.lazy import pylab
from sequana.lazy import pandas as pd
from sequana.lazy import numpy as np
from sequana.viz.core import VizInput2D

import colorlog
logger = colorlog.getLogger(__name__)




__all__ = ["Hist2D"]


class Hist2D(VizInput2D):
    """2D histogram

    .. plot::
        :include-source:
        :width: 80%

        from numpy import random
        from sequana.viz import hist2d
        X = random.randn(10000)
        Y = random.randn(10000)
        h = hist2d.Hist2D(X,Y)
        h.plot(bins=100, contour=True)

    """
    def __init__(self, x, y=None, verbose=False):
        """.. rubric:: constructor

        :param x: an array for X values. See :class:`~sequana.viz.core.VizInput2D` for details.
        :param y: an array for Y values. See :class:`~sequana.viz.core.VizInput2D` for details.

        """
        super(Hist2D, self).__init__(x=x, y=y, verbose=verbose)

    def plot(self, bins=100, cmap="hot_r", fontsize=10, Nlevels=4,
        xlabel=None, ylabel=None, norm=None, range=None, normed=False,
        colorbar=True, contour=True, grid=True, **kargs):
        """plots histogram of mean across replicates versus coefficient variation

        :param int bins: binning for the 2D histogram (either a float or list
            of 2 binning values).
        :param cmap: a valid colormap (defaults to hot_r)
        :param fontsize: fontsize for the labels
        :param int Nlevels: must be more than 2
        :param str xlabel: set the xlabel (overwrites content of the dataframe)
        :param str ylabel: set the ylabel (overwrites content of the dataframe)
        :param norm: set to 'log' to show the log10 of the values.
        :param normed: normalise the data
        :param range: as in pylab.Hist2D : a 2x2 shape [[-3,3],[-4,4]]
        :param contour: show some contours (default to True)
        :param bool grid: Show unerlying grid (defaults to True)

        If the input is a dataframe, the xlabel and ylabel will be populated
        with the column names of the dataframe.

        """
        X = self.df[self.df.columns[0]].values
        Y = self.df[self.df.columns[1]].values
        if len(X) > 10000:
            logger.info("Computing 2D histogram. Please wait")

        pylab.clf()
        if norm == 'log':
            from matplotlib import colors
            res = pylab.hist2d(X, Y, bins=bins, density=normed,
               cmap=cmap, norm=colors.LogNorm())
        else:
            res = pylab.hist2d(X, Y, bins=bins, cmap=cmap,
                    density=normed, range=range)

        if colorbar is True:
            pylab.colorbar()

        if contour:
            try:
                bins1 = bins[0]
                bins2 = bins[1]
            except:
                bins1 = bins
                bins2 = bins

            X, Y = pylab.meshgrid(res[1][0:bins1], res[2][0:bins2])
            if contour:
                if res[0].max().max() < 10 and norm == 'log':
                    pylab.contour(X, Y, res[0].transpose())
                else:
                    levels = [round(x) for x in
                            pylab.logspace(0, pylab.log10(res[0].max().max()), Nlevels)]
                    pylab.contour(X, Y, res[0].transpose(), levels[2:])
                #pylab.clabel(C, fontsize=fontsize, inline=1)

        if ylabel is None:
            ylabel = self.df.columns[1]
        if xlabel is None:
            xlabel = self.df.columns[0]

        pylab.xlabel(xlabel, fontsize=fontsize)
        pylab.ylabel(ylabel, fontsize=fontsize)

        if grid is True:
            pylab.grid(True)

        return res


