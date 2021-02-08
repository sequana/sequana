# -*- coding: utf-8 -*-
#
#  This file is part of Sequana software
#
#  Copyright (c) 2016 - Sequana Development Team
#
#  File author(s):
#      Thomas Cokelaer <thomas.cokelaer@pasteur.fr>
#
#  Distributed under the terms of the 3-clause BSD license.
#  The full license is in the LICENSE file, distributed with this software.
#
#  website: https://github.com/sequana/sequana
#  documentation: http://sequana.readthedocs.io
#
##############################################################################
"""Imshow utility"""
from sequana.viz.core import VizInputSquare
from sequana.lazy import pylab
from sequana.lazy import pandas as pd


__all__ = ['Imshow']


class Imshow(VizInputSquare):
    """Wrapper around the matplotlib.imshow function

    Very similar to matplotlib but set interpolation to None, and aspect
    to automatic and accepts input as a dataframe, in whic case
    x and y labels are set automatically.


    .. plot::
        :width: 80%
        :include-source:

        import pandas as pd
        data = dict([ (letter,np.random.randn(10)) for letter in 'ABCDEFGHIJK'])
        df = pd.DataFrame(data)

        from sequana.viz import Imshow
        im = Imshow(df)
        im.plot()

    """
    def __init__(self, x, verbose=True):
        """.. rubric:: constructor

        :param x: input dataframe (or numpy matrix/array). 
            Must be squared.

        """
        super(Imshow, self).__init__(x, verbose=verbose)

    def plot(self, interpolation='None', aspect='auto', cmap='hot', tight_layout=True,
        colorbar=True, fontsize_x=None, fontsize_y=None, rotation_x=90,
        xticks_on=True, yticks_on=True, **kargs):
        """wrapper around imshow to plot a dataframe

        :param interpolation: set to None
        :param aspect: set to 'auto'
        :param cmap: colormap to be used.
        :param tight_layout:
        :param colorbar: add a colobar (default to True)
        :param fontsize_x: fontsize on xlabels
        :param fontsize_y: fontsize on ylabels
        :param rotation_x: rotate labels on xaxis
        :param xticks_on: switch off the xticks and labels
        :param yticks_on: switch off the yticks and labels

        """

        data = self.df
        pylab.clf()
        pylab.imshow(data, interpolation=interpolation, aspect=aspect, cmap=cmap, **kargs)

        if fontsize_x == None:
            fontsize_x = 16 #FIXME use default values
        if fontsize_y == None:
            fontsize_y = 16 #FIXME use default values

        if yticks_on is True:
            pylab.yticks(range(0, len(data.index)), data.index, 
                fontsize=fontsize_y)
        else:
            pylab.yticks([])
        if xticks_on is True:
            pylab.xticks(range(0, len(data.columns[:])), data.columns, 
                fontsize=fontsize_x, rotation=rotation_x)
        else:
            pylab.xticks([])

        if colorbar is True:
            pylab.colorbar()

        if tight_layout:
            pylab.tight_layout()

        # For some reasons, in newest version of python/mpl, this is required
        # for ylim, not for xlim
        y1,y2 = pylab.ylim()
        pylab.ylim([y1+0.5, y2-0.5])


