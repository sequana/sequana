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


from sequana.viz import Imshow
from sequana.lazy import numpy as np
from sequana.lazy import pandas as pd


__all__ = ['ANOVA']


class ANOVA(object):
    """DRAFT

    Testing if 3(+) population means are all equal.

    Looks like the group are different, visually, and naively. 

    .. plot:: 
        :include-source:

        from pylab import *
        from sequana.viz import ANOVA
        import pandas as pd

        A = normal(0.5,size=10000)
        B = normal(0.25, size=10000)
        C = normal(0, 0.5,size=10000)
        df = pd.DataFrame({"A":A, "B":B, "C":C})
        a = ANOVA(df)
        print(a.anova())
        a.imshow_anova_pairs()

    """

    def __init__(self, df):

        self.df = df.copy()

    def anova(self):
        """Perform one-way ANOVA.

        The one-way ANOVA tests the null hypothesis that two or more groups have
        the same population mean.  The test is applied to samples from two or
        more groups. Since we ar using a dataframe, vector length are identical.

        return: the F value (test itself), and its p-value
        """
        import scipy.stats
        F, P = scipy.stats.f_oneway(*[self.df[x] for x in self.df.columns])
        return F, P

    def imshow_anova_pairs(self, log=True, **kargs):
        import scipy.stats
        N = len(self.df.columns)

        # could use a dataframe straight way ?
        res = np.ones((N, N))
        for i,col1 in enumerate(self.df.columns):
            for j,col2 in enumerate(self.df.columns):
                d1 = self.df[col1]
                d2 = self.df[col2]
                F, P = scipy.stats.f_oneway(*[d1, d2])
                res[i][j] = P
        df = pd.DataFrame(res, index=self.df.columns, columns=self.df.columns)
        #FIXME: may have na, which are set to 1
        df = df.fillna(1)
        if log is True:
            Imshow(-np.log10(df)).plot(**kargs)
        else:
            Imshow(df).plot(**kargs)
        return df

