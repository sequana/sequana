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


__all__ = ["Cluster"]


class Cluster:
    """

    Input must be a matrix in the form of a pandas DataFrame. Each column is a
    sample. Sample names are the columns' names. colors are set to red for all
    samples but user can provide a mapping of columns' names and a color.


    ::

        c = Cluster(data, colors={"A": "r", "B": "g"}


    """

    def __init__(self, data, colors={}):
        """.. rubric:: constructor

        :param data: a dataframe; Each column being a sample.
        :param colors: a mapping of column/sample name a color

        """
        self.df = data
        self.labels = data.columns
        self.colors = {x: "r" for x in self.labels}
        for k, v in colors.items():
            self.colors[k] = v

        from sklearn.preprocessing import StandardScaler

        self.scaler = StandardScaler()

    def scale_data(self, transform_method="log", max_features=500):
        """

        - Replace zeros with 1 (avoid log issue)
        - transform the data using log10 or anscombe transform
        - scale the data using the scaler attribute (standard scaler by default)

        """
        assert transform_method in ["log", "anscombe"]
        # normalise the data

        # First, we transform the data
        data = self.df.copy()
        data = data.replace(0, 1)
        self.data = data
        if transform_method == "log":
            data = pylab.log10(data)
        elif transform_method == "anscombe":
            from sequana.vst import VST

            data = VST.anscombe(data)

        # then we keep only the first N most dispersed features
        tokeep = data.std(axis=1).sort_values(ascending=False).index[0:max_features]
        data = data.loc[tokeep]
        data = self.scaler.fit_transform(data)
        return data, tokeep

    def _plot(
        self, Xr, pca=None, pc1=0, pc2=1, colors=None, show_labels=True, fontsize=10
    ):

        if colors is None:
            colors = [self.colors[k] for k in self.labels]
            if len(colors) != len(Xr):
                colors = ["r"] * len(Xr[:, 0])
        else:
            for k in self.labels:
                if k not in colors.keys():
                    logger.warning(
                        "No key color for this sample: {}. Set to red".format(k)
                    )
                    colors[k] = "r"
            colors = [colors[k] for k in self.labels]

        pylab.scatter(Xr[:, pc1], Xr[:, pc2], c=colors)
        ax = pylab.gca()
        X1, X2 = pylab.xlim()
        dX = X2 - X1
        pylab.xlim([X1 + X1 * 0.05, X2 + X2 * 0.05])

        Y1, Y2 = pylab.ylim()
        dY = Y2 - Y1
        pylab.ylim([Y1 + Y1 * 0.05, Y2 + Y2 * 0.05])

        count = 0
        if fontsize == 0:
            show_labels = 0
        if show_labels:
            for x, y in zip(Xr[:, pc1], Xr[:, pc2]):
                x += dX / 40
                y += dY / 40
                ax.annotate(
                    self.labels[count], (x, y), color=colors[count], fontsize=fontsize
                )
                count += 1
                if count > 100:
                    break
        if pca:
            pylab.xlabel(
                "PC{} ({}%)".format(
                    pc1 + 1, round(pca.explained_variance_ratio_[pc1] * 100, 2)
                ),
                fontsize=12,
            )
            pylab.ylabel(
                "PC{} ({}%)".format(
                    pc2 + 1, round(pca.explained_variance_ratio_[pc2] * 100, 2)
                ),
                fontsize=12,
            )
        pylab.grid(True)
