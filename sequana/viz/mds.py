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

from sequana.lazy import pylab

import colorlog
logger = colorlog.getLogger(__name__)

from sequana.viz import clusterisation
__all__ = ['Isomap']


class MDS(clusterisation.Cluster):
    """

    .. plot::
        :include-source:

        from sequana.viz.mds import MDS
        from sequana import sequana_data
        import pandas as pd

        data = sequana_data("test_pca.csv")
        df = pd.read_csv(data)
        df = df.set_index("Id")
        p = PCA(df, colors={
            "A1": 'r', "A2": 'r', 'A3': 'r',
            "B1": 'b', "B2": 'b', 'B3': 'b'})
        p.plot(n_components=2)

    """
    def __init__(self, data, colors={}):
        super(MDS, self).__init__(data, colors)

    def plot(self, n_components=2, transform="log", switch_x=False,
            switch_y=False, switch_z=False, colors=None,
            max_features=500, show_plot=True, clf=True, show_labels=True):
        """

        :param n_components: at number starting at 2 or a value below 1
            e.g. 0.95 means select automatically the number of components to
            capture 95% of the variance
        :param transform: can be 'log' or 'anscombe', log is just log10. count
            with zeros, are set to 1
        """
        from sklearn.manifold import MDS
        import numpy as np

        if clf: pylab.clf()

        data, kept = self.scale_data(transform_method=transform, max_features=max_features)

        mds = MDS(n_components=n_components)
        Xr = mds.fit_transform(data.T)
        self.Xr = Xr

        if switch_x:
            Xr[:,0] *= -1
        if switch_y:
            Xr[:,1] *= -1
        if switch_z:
            Xr[:,2] *= -1

        # PC1 vs PC2
        if show_plot:
            pylab.figure(1)
            self._plot(Xr, pca=None, pc1=0,pc2=1, colors=colors,
                    show_labels=show_labels)

        if n_components >=3:
            if show_plot:
                pylab.figure(2)
                self._plot(Xr, pca=None, pc1=0,pc2=2, colors=colors,
                    show_labels=show_labels)
                pylab.figure(3)
                self._plot(Xr, pca=None, pc1=1,pc2=2, colors=colors,
                    show_labels=show_labels)
        return mds



