#
#  This file is part of Sequana software
#
#  Copyright (c) 2016-2021 - Sequana Development Team
#
#
#  Distributed under the terms of the 3-clause BSD license.
#  The full license is in the LICENSE file, distributed with this software.
#
#  website: https://github.com/sequana/sequana
#  documentation: http://sequana.readthedocs.io
#
##############################################################################

import colorlog
from adjustText import adjust_text

from sequana.lazy import numpy as np
from sequana.lazy import pylab
from sequana.viz import clusterisation

logger = colorlog.getLogger(__name__)


__all__ = ["PCA"]


class PCA(clusterisation.Cluster):
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


    From R, a PCA is selecting the first 500 features based on variance.


    """

    def __init__(self, data, colors={}):
        super(PCA, self).__init__(data, colors)

    def plot_pca_vs_max_features(self, step=100, n_components=2, transform="log"):
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
        assert n_components in [2, 3, 4]
        N = len(self.df)
        if step > N:
            step = N

        # We start with at least 5 features
        X = range(10, N, step)
        Y = []
        for i, x in enumerate(X):
            res = self.plot(n_components=n_components, max_features=x, show_plot=False, transform=transform)
            Y.append(res)

        sub = n_components
        pylab.subplot(sub, 1, 1)
        pylab.plot(X, [y[0] * 100 for y in Y])
        pylab.ylabel("PC1 (%)")
        pylab.subplot(sub, 1, 2)
        pylab.plot(X, [y[1] * 100 for y in Y])
        pylab.ylabel("PC2 (%)")
        if sub >= 3:
            pylab.subplot(sub, 1, 3)
            pylab.plot(X, [y[2] * 100 for y in Y])
            pylab.ylabel("PC3 (%)")
        if sub >= 4:
            pylab.subplot(sub, 1, 4)
            pylab.plot(X, [y[3] * 100 for y in Y])
            pylab.ylabel("PC4 (%)")

    def plot(
        self,
        n_components=2,
        transform="log",
        switch_x=False,
        switch_y=False,
        switch_z=False,
        colors=None,
        max_features=500,
        show_plot=True,
        fontsize=10,
        adjust=True,
    ):
        """

        :param n_components: at number starting at 2 or a value below 1
            e.g. 0.95 means select automatically the number of components to
            capture 95% of the variance
        :param transform: can be 'log' or 'anscombe', log is just log10. count
            with zeros, are set to 1
        """
        from sklearn.decomposition import PCA

        pylab.clf()

        pca = PCA(n_components)

        # scale the data
        data = self.scale_data(transform_method=transform)

        # keep only top variable features
        tokeep = data.std(axis=1).sort_values(ascending=False).index[0:max_features]
        data = data.loc[tokeep]

        pca.fit(data.T)

        if transform == "standard":
            # If scale data used a scaler, then we need to use it for the transformation
            Xr = pca.transform(self.scaler.fit_transform(self.df.loc[tokeep].T))
        else:
            # otherwise, noting to do
            Xr = pca.transform(self.df.loc[tokeep].T)
        self.Xr = Xr

        if switch_x:
            Xr[:, 0] *= -1
        if switch_y:
            Xr[:, 1] *= -1
        if switch_z:
            Xr[:, 2] *= -1

        def adjust_func():
            texts = [x for x in pylab.gca().get_children() if "text.Annotation" in str(x.__class__)]
            adjust_text(texts)

        # PC1 vs PC2
        if show_plot:
            pylab.figure(1)
            self._plot(Xr, pca=pca, pc1=0, pc2=1, colors=colors, fontsize=fontsize)
            if adjust:
                adjust_func()

        if len(pca.explained_variance_ratio_) >= 3:
            if show_plot:
                pylab.figure(2)
                self._plot(Xr, pca=pca, pc1=0, pc2=2, colors=colors, fontsize=fontsize)
                if adjust:
                    adjust_func()

                pylab.figure(3)
                self._plot(Xr, pca=pca, pc1=1, pc2=2, colors=colors, fontsize=fontsize)
                if adjust:
                    adjust_func()

        return pca.explained_variance_ratio_
