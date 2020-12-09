# -*- coding: utf-8 -*-
#
#  This file is part of Sequana software
#
#  Copyright (c) 2016-2020 - Sequana Development Team
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
"""Heatmap and dendograms"""
import matplotlib
import pylab
import scipy.cluster.hierarchy as hierarchy
import scipy.spatial.distance as distance
import numpy as np # get rid of this dependence

import easydev
import colormap
from sequana.viz.linkage import Linkage

__all__ = ['Dendogram']


class Dendogram(Linkage):
    """dendograms of an input matrix

    .. plot::
        :include-source:
        :width: 80%

        from sequana.viz import heatmap, dendogram
        df = heatmap.get_heatmap_df()
        h = dendogram.Dendogram(df)
        h.plot()


    You should scale the data before::

        from sequana.viz.clusterisation import Clusterisation
        scaled, index = Clusterisation(data).scale_data()
        import pandas as pd
        df = pd.DataFrame(scaled)
        df.index = index
        df.columns = data.columns

        g = Dendogram(df.T)
        g.plot()

    """
    def __init__(self, data=None, method='complete', 
                 metric='euclidean',
                 cmap='yellow_black_blue',
                 col_side_colors=None, side_colors=None,
                 verbose=True, horizontal=True
                 ):
        """.. rubric:: constructor

        :param data: a dataframe or possibly a numpy matrix.
        :param method: complete by default
        :param metric: euclidean by default
        :param cmap: colormap. any matplotlib accepted or combo of colors as
            defined in colormap package (pypi)
        :param col_side_colors:
        :param side_colors:

        """
        # should be a copy since it may be reshuffled ?
        try:
            if data is None and verbose is True:
                print("No data provided, please fill the `df` attribute manually")
            elif data is None:
                pass
            else:
                self._df = data.copy()
        except AttributeError as err:
            print("input must be a pandas data frame or numpy matrix")
            raise(err)

        self._method = method
        self._metric = metric
        self.horizontal = True

        # some default parameters
        self.cluster_criterion = 'distance'
        self.params = easydev.AttrDict()
        self.params.side_colors = ['r', 'g', 'b', 'y', 'w', 'k', 'm']
        self.params.cmap = cmap

        self.category = {}

        if side_colors:
            self.params.side_colors = side_colors

    def _get_df(self):
        return self._df
    def _set_df(self, data):
        self._df = data.copy()
    df = property(_get_df, _set_df)
    frame = property(_get_df, _set_df)

    def _get_method(self):
        return self._method
    def _set_method(self, value):
        self.check_method(value)
        self._method = value
    method = property(_get_method, _set_method)

    def _get_metric(self):
        return self._metric
    def _set_metric(self, value):
        self.check_metric(value)
        self._metric = value
    metric = property(_get_metric, _set_metric)

    def plot(self, num=1, cmap=None, colorbar=True, 
             figsize=(12, 8),
             fontsize=None
             ):
        """

        Using as input::

            df = pd.DataFrame({'A':[1,0,1,1],
                               'B':[.9,0.1,.6,1],
                            'C':[.5,.2,0,1],
                            'D':[.5,.2,0,1]})

        .. plot::
            :include-source:
            :width: 80%

            from sequana.viz import heatmap
            df = heatmap.get_heatmap_df()
            h = heatmap.Heatmap(df)
            h.category_row['A'] = 1
            h.category_row['C'] = 1
            h.category_row['D'] = 2
            h.category_row['B'] = 2
            h.plot()


        """
        # save all parameters in a dict
        layout = {}

        if cmap is None:
            cmap = self.params.cmap
        try:cmap = colormap.cmap_builder(cmap)
        except:pass

        # keep track of row and column names for later.
        header = self.frame.index

        # FIXME something clever for the fontsize
        if len(header) > 100 or len(header) > 100:
            matplotlib.rcParams['font.size'] = 6
        if len(header) > 50 or len(header) > 50:
            matplotlib.rcParams['font.size'] = 7
        if len(header) > 30 or len(header) > 30:
            matplotlib.rcParams['font.size'] = 8
        else:
            matplotlib.rcParams['font.size'] = 12
        if fontsize:
            matplotlib.rcParams['font.size'] = fontsize

        # scaling min/max range

        # Scale the figure window size #
        fig = pylab.figure(num=num, figsize=figsize)
        fig.clf()

        Y = self.linkage(self.frame, self.method, self.metric )

        Z = hierarchy.dendrogram(Y, orientation='right',
             color_threshold=0,
            above_threshold_color="k", distance_sort="descending")
        ind1 = hierarchy.fcluster(Y, 0.7 * max(Y[:,2]), self.cluster_criterion)

        # apply the clustering for the array-dendrograms to the actual matrix data
        idx1 = Z['leaves']

        # Rearrange the data frame in the order of the dendogram
        self.frame = self.frame.iloc[idx1,:]
        ticks = pylab.yticks()[0]
        pylab.yticks(ticks, self.frame.index)
        pylab.tight_layout()

        # reorder the flat cluster to match the order of the leaves the dendrogram
        ind1 = ind1[idx1]

        if self.category:
            gca = pylab.gca()
            X, Y = gca.get_position().get_points()
            f = pylab.gcf()
            ax = f.add_axes([X[0], X[1], 0.02, Y[1]-X[1]])

            category= [self.category[x] for x in self.df.index]
            dr = np.array(category, dtype=int)
            dr.shape = (len(category),1)
            cmap_r = matplotlib.colors.ListedColormap(self.params.side_colors)
            ax.matshow(dr, aspect='auto', origin='lower', cmap=cmap_r)
            ax.set_xticks([])
            ax.set_yticks([])
        


