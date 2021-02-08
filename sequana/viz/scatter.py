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

""".. rubric:: Scatter plots

:author: Thomas Cokelaer
"""
from sequana.lazy import pylab
from sequana.lazy import pandas as pd
from sequana.viz.core import VizInput2D

__all__ = ["ScatterHist"]



class ScatterHist(VizInput2D):
    """Scatter plots and histograms


    """
    def __init__(self, x, y=None, verbose=True):
        """.. rubric:: constructor

        :param x: if x is provided, it should be a dataframe with 2 columns. The
            first one will be used as your X data, and the second one as
            the Y data
        :param y:
        :param verbose:


        """
        super(ScatterHist, self).__init__(x,y,verbose)

    def plot(self,
        kargs_scatter={'s':20, 'c':'b'},
        kargs_grids={},
        kargs_histx={},
        kargs_histy={},
        scatter_position='bottom left',
        width=.5,
        height=.5,
        offset_x=.10,
        offset_y=.10,
        gap=0.06,
        facecolor='lightgrey',
        grid=True,
        show_labels=True,
        **kargs):
        """Scatter plot of set of 2 vectors and their histograms.

        :param x: a dataframe or a numpy matrix (2 vectors) or a list of 2 items,
            which can be a mix of list or numpy array.
            if **size** and/or **color** are found in the columns dataframe,
            those columns will be used in the scatter plot. kargs_scatter keys **c**
            and **s** will then be ignored. If a list of lists, **x** will be the first row
            and **y** the second row.
        :param y: if x is a list or an array, then y must also be provided as
            a list or an array
        :param kargs_scatter: a dictionary with pairs of key/value accepted by
            matplotlib.scatter function. Examples is a list of colors or a list
            of sizes as shown in the examples below.
        :param kargs_grid: a dictionary with pairs of key/value accepted by
            the maplotlib.grid (applied on histogram and axis at the same time)
        :param kargs_histx: a dictionary with pairs of key/value accepted by the
            matplotlib.histogram
        :param kargs_histy: a dictionary with pairs of key/value accepted by the
            matplotlib.histogram
        :param kargs: other optional parameters are **hold**, **facecolor**.
        :param scatter_position: can be 'bottom right/bottom left/top left/top right'
        :param width: width of the scatter plot (value between 0 and 1)
        :param height: height of the scatter plot (value between 0 and 1)
        :param offset_x:
        :param offset_y:
        :param gap: gap between the scatter and histogram plots.
        :param grid: defaults to True

        :return: the scatter, histogram1 and histogram2 axes.

        .. plot::
            :include-source:
            :width: 80%

            import pylab
            import pandas as pd
            X = pylab.randn(1000)
            Y = pylab.randn(1000)
            df = pd.DataFrame({'X':X, 'Y':Y})

            from sequana.viz import ScatterHist
            ScatterHist(df).plot()


        .. plot::
            :include-source:
            :width: 80%

            from sequana.viz import ScatterHist
            ScatterHist(x=[1,2,3,4], y=[3,5,6,4]).plot(
                kargs_scatter={
                    's':[200,400,600,800],
                    'c': ['red', 'green', 'blue', 'yellow'],
                    'alpha':0.5},
                kargs_histx={'color': 'red'},
                kargs_histy={'color': 'green'})


        .. seealso:: `notebook <http://nbviewer.ipython.org/github/sequana/sequana/blob/master/notebooks/viz/scatter.ipynb>`__
        """
        df = self.df
        try:
            kargs_scatter['s'] = df['size']
        except:
            pass
        try:
            kargs_scatter['c'] = df['color']
        except:
            pass


        if kargs.get("hold", False) is False:
            pylab.clf()

        W = width
        H = height
        if scatter_position == 'bottom left':
            X0 = offset_x
            Y0 = offset_y
            Xoff = X0 + W + gap
            Yoff = Y0 + H + gap
            Wh = 1 - offset_x*2 - W - gap
            Hh = 1 - offset_y*2 - H - gap
        elif scatter_position == 'bottom right':
            Wh = 1 - offset_x*2 - W - gap
            Hh = 1 - offset_y*2 - H - gap
            X0 = offset_x + Wh +gap
            Y0 = offset_y
            Xoff = offset_x
            Yoff = Y0 + H + gap
        elif scatter_position == 'top right':
            Wh = 1 - offset_x*2 - W - gap
            Hh = 1 - offset_y*2 - H - gap
            X0 = offset_x + Wh +gap
            Y0 = offset_y + Hh + gap
            Xoff = offset_x
            Yoff = offset_y
        elif scatter_position == 'top left':
            Wh = 1 - offset_x*2 - W - gap
            Hh = 1 - offset_y*2 - H - gap
            X0 = offset_x
            Y0 = offset_y + Hh + gap
            Xoff = offset_x + W + gap
            Yoff = offset_y #Y0 #+ H + gap
        else: #pragma: no cover
            raise ValueError("scatter_position must be 'top left', 'top right', 'bottom left', 'bottom right'")

        facecolor = kargs.get('facecolor', 'lightgrey')



        ax_scatter = pylab.axes((X0, Y0, W, H), facecolor=facecolor, xscale='linear',
            yscale='linear')#, xticks='auto', yticks='auto')

        if show_labels:
            ax_scatter.set_xlabel(self.xy_names[0])
            ax_scatter.set_ylabel(self.xy_names[1])
        ax_hist_x = pylab.axes((X0, Yoff, W, Hh), facecolor=facecolor, xscale='linear',
            yscale='linear')#, xticks='auto', yticks='auto')
        ax_hist_y = pylab.axes((Xoff, Y0, Wh, H), facecolor=facecolor, xscale='linear',
            yscale='linear')#, xticks='auto', yticks='auto')

        # move ticks on axis  if needed
        ax_hist_x.xaxis.set_ticks_position('top')
        if scatter_position == 'bottom left':
            ax_scatter.yaxis.set_ticks_position('left')
            ax_hist_x.yaxis.set_ticks_position('right')
        elif scatter_position == 'bottom right':
            ax_hist_y.yaxis.set_ticks_position('left')
        elif scatter_position == 'top right':
            ax_scatter.xaxis.set_ticks_position('top')
            ax_scatter.yaxis.set_ticks_position('right')
            ax_hist_y.yaxis.set_ticks_position('left')
            ax_hist_x.xaxis.set_ticks_position('bottom')
        elif scatter_position == 'top left':
            ax_scatter.xaxis.set_ticks_position('top')
            ax_hist_y.yaxis.set_ticks_position('right')
            ax_hist_x.xaxis.set_ticks_position('bottom')
        else: #pragma: no cover
            raise ValueError("scatter_position must be 'top left', 'top right', 'bottom left', 'bottom right'")

        ax_scatter.scatter(df.x, df.y, **kargs_scatter)
        ax_hist_x.hist(df.x, **kargs_histx)
        # fixme: user may not want that ?
        kargs_histy['orientation'] = 'horizontal'
        ax_hist_y.hist(df.y, **kargs_histy)
        # I tried c.set_xticks but rotation could not be found
        pylab.xticks(ax_hist_y.get_xticks(), rotation=90)

        # grid
        if grid is True:
            ax_scatter.grid(b=grid, which='major', axis='both', **kargs_grids)
            ax_hist_x.grid(b=grid, which='major', axis='both', **kargs_grids)
            ax_hist_y.grid(b=grid, which='major', axis='both', **kargs_grids)

        return (ax_scatter, ax_hist_x, ax_hist_y)





