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
from sequana import logger

__all__ = ['PCA']


class PCA():
    """

    .. plot::
        :include-source:

        from sequana.viz.pca import PCA
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

        self.df = data
        self.labels = data.columns
        self.colors = {x:'r' for x in self.labels}
        for k,v in colors.items():
            self.colors[k] = v

    def plot_pca_vs_max_features(self, max_features=100, n_components=2):
        """

        .. plot::
            :include-source:

            from sequana.viz.pca import PCA
            from sequana import sequana_data
            import pandas as pd

            data = sequana_data("test_pca.csv")
            df = pd.read_csv(data)
            df = df.set_index("Id")

            p = PCA(df)
            p.plot_pca_vs_max_features()

        """
        assert n_components in [2,3]
        N = len(self.df)
        if max_features > N:
            max_features = N
        # We start with at least 5 features
        X = range(5, N, max_features)
        Y = [self.plot(n_components=n_components, max_features=x, show_plot=False) for x in X]
        sub = n_components
        pylab.subplot(sub,1,1)
        pylab.plot(X, [y[0]*100 for y in Y])
        pylab.ylabel("PC1 (%)")
        pylab.subplot(sub,1,2)
        pylab.plot(X, [y[1]*100 for y in Y])
        pylab.ylabel("PC2 (%)")
        if sub == 3:
            pylab.subplot(sub,1,3)
            pylab.plot(X, [y[2]*100 for y in Y])
            pylab.ylabel("PC3 (%)")

    def plot(self, n_components=2, transform="log", switch_x=False,
            switch_y=False, switch_z=False, colors=None,
            max_features=500, show_plot=True):
        """

        :param n_components: at number starting at 2 or a value below 1
            e.g. 0.95 means select automatically the number of components to
            capture 95% of the variance
        :param transform: can be 'log' or 'anscombe', log is just log10. count
            with zeros, are set to 1
        """
        assert transform in ['log', 'anscombe']
        
        from sklearn.preprocessing import StandardScaler
        from sklearn.decomposition import PCA
        import numpy as np

        pylab.clf()
        pca = PCA(n_components)

        # normalise the data
        scaler = StandardScaler()

        # First, we transform the data
        data = self.df.copy()
        data = data.replace(0, 1)
        self.data = data
        if transform == "log":
            data = pylab.log10(data)
        elif transform == "anscombe":
            from sequana.vst import VST
            data = VST.anscombe(data)

        # then we keep only the first N most dispersed features
        tokeep = data.std(axis=1).sort_values(ascending=False).index[0:max_features]
        data = data.loc[tokeep]

        data = scaler.fit_transform(data)

        pca.fit(data.T)

        Xr = pca.transform(scaler.fit_transform(self.df.loc[tokeep].T))
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
            self._plot(Xr, pca, 0,1, colors=colors)

        if len(pca.explained_variance_ratio_) >= 3:
            if show_plot:
                pylab.figure(2)
                self._plot(Xr, pca, 0,2, colors=colors)
                pylab.figure(3)
                self._plot(Xr, pca, 1,2, colors=colors)

        return pca.explained_variance_ratio_

    def _plot(self, Xr, pca, pc1=0, pc2=1, colors=None):
        if colors is None:
            colors = [self.colors[k] for k in self.labels]
            if len(colors) != len(Xr):
                colors = ["r"] * len(Xr[:,0])
        else:
            for k in self.labels:
                if k not in colors.keys():
                    logger.warning("No key color for this sample: {}. Set to red".format(k))
                    colors[k] = "r"
            colors = [colors[k] for k in self.labels]

        pylab.scatter(Xr[:,pc1], Xr[:,pc2], c=colors)
        ax = pylab.gca()
        X1, X2 = pylab.xlim()
        dX = X2 - X1
        pylab.xlim([X1 + X1*0.05, X2 + X2*0.05])

        Y1, Y2 = pylab.ylim()
        dY = Y2 - Y1
        pylab.ylim([Y1 + Y1*0.05, Y2 + Y2*0.05])

        count = 0
        for x,y in zip(Xr[:,pc1], Xr[:,pc2]):
            x += dX / 40
            y += dY / 40
            ax.annotate(self.labels[count], (x,y))
            count += 1
            if count > 100: 
                break
        pylab.xlabel("PC{} ({}%)".format(pc1+1,
            round(pca.explained_variance_ratio_[pc1]*100, 2)))
        pylab.ylabel("PC{} ({}%)".format(pc2+1,
            round(pca.explained_variance_ratio_[pc2]*100, 2)))
        pylab.grid(True)


