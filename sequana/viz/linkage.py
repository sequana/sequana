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

"""Heatmap and dendograms"""
import warnings

import scipy.cluster.hierarchy as hierarchy
import scipy.spatial.distance as distance

__all__ = ["Linkage"]


class Linkage(object):
    """Linkage used in other tools such as Heatmap"""

    methods = ["single", "complete", "average", "weighted", "centroid", "median", "ward"]
    metrics = [
        "braycurtis",
        "canberra",
        "chebyshev",
        "cityblock",
        "correlation",
        "cosine",
        "dice",
        "euclidean",
        "hamming",
        "jaccard",
        "jensenshannon",
        "kulsinski",
        "mahalanobis",
        "matching",
        "minkowski",
        "rogerstanimoto",
        "russellrao",
        "seuclidean",
        "sokalmichener",
        "sokalsneath",
        "sqeuclidean",
        "yule",
    ]

    def __init__(self):
        """.. rubric:: constructor

        :param data: a dataframe or possibly a numpy matrix.

        """
        pass

    def check_metric(self, value):
        if value not in self.metrics:
            raise ValueError(f"value {value} must be in {self.metrics}")

    def check_method(self, value):
        # None is possible
        # in R, in addition to single, complete, average, centroid,
        # median and ward
        # there are  ward.D, wardD2 and mcquitty
        # default is complete
        if value not in self.methods:
            raise ValueError(f"value {value} must be in {self.methods}")

    def linkage(self, df, method, metric):
        self.check_metric(metric)
        self.check_method(method)
        d = distance.pdist(df)
        D = distance.squareform(d)

        # hierarchy.ClusterWarning
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            Y = hierarchy.linkage(D, method=method, metric=metric)
            return Y
