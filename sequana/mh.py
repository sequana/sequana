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



from random import uniform, gauss
from sequana.lazy import numpy as np
from sequana.lazy import pylab

import colorlog
logger = colorlog.getLogger(__name__)




__all__ = ["MetropolisHasting"]


class MetropolisHasting():
    """

    .. plot::

        from sequana.mh import MetropolisHasting
        m = MetropolisHasting()
        m.Xtarget = [0.   ,  0.005,  0.01 ,  0.016,  0.021,  0.027,  0.032,  0.037,
            0.043,  0.048,  0.054,  0.059,  0.065,  0.07 ,  0.075,  0.081,
            0.086,  0.092,  0.097,  0.103]
        m.Ytarget = [83, 315,  611, 675, 1497, 5099, 7492, 2797, 842, 334, 
            117, 63, 33, 22, 11, 3, 3, 1,  0,  2]
        vec = m.simulate(100000)
        m.check(bins=100)


    .. warning: be aware of border effects. For instance if the profile does
        not go to zero at the lower bound or upper_bound, then the final
        histogram may be biased on the boundaries. One would need to incrase the
        boundaries by a range larger than the step used in the proposal / jump
        function and remove the data outside of the expected boundaries
        manually. 

    """
    def __init__(self):
        self.burning = 20000
        self.lower_bound = 0000
        self.upper_bound = 100000
        self.aprob = []
        self._Xtarget = None

    def _set_x(self, X):
        self._Xtarget = np.array(X)
        self.lower_bound = min(self._Xtarget)
        self.upper_bound = max(self._Xtarget)
    def _get_x(self):
        return self._Xtarget
    Xtarget = property(_get_x, _set_x)

    def _set_y(self, Y):
        self._Ytarget = np.array(Y)
    def _get_y(self):
        return self._Ytarget
    Ytarget = property(_get_y, _set_y)

    def simulate(self, n=100000, burning=20000, step=None, x0=None):
        if step is None:
            self.step = (self.upper_bound - self.lower_bound) / 100.
            step = self.step

        if x0 is None:
            self.x0 = (self.upper_bound - self.lower_bound) / 2 + self.lower_bound

        self.aprob = []

        # function target profile
        NS = self.Ytarget / sum(self.Ytarget)
        sdnorm = lambda x: np.interp(x, self.Xtarget, NS)

        x = self.x0

        vec = [x] # starting seed

        # a gaussian jump centered on 0 is used as a random inovation 
        # if the candidate is outside of the boundaries, we try another 
        # candidate
        def jumper(x):
            #jump = uniform(-step, step)
            jump = gauss(0, step)
            xprime = x + jump
            while xprime < self.lower_bound or xprime>self.upper_bound:
                jump = gauss(0, step)
                xprime = x + jump
            return xprime

        for i in range(1, n*2+burning):
            xprime = jumper(x)

            aprob = min([1., sdnorm(xprime)/sdnorm(x)]) #acceptance probability
            u = uniform(0, 1)
            if u < aprob:
                x = xprime
                vec.append(x)
            self.aprob.append(aprob)

            if len(vec) == n + burning:
                break
        self.burning_vector = vec[0:burning]
        self.vec = vec[burning:]

        return vec[burning:]

    def diagnostics(self, bins=60, clear=True):
        if clear: pylab.clf()

        pylab.subplot(3,1,1)
        pylab.hist(self.aprob, bins=bins)
        pylab.title("Acceptation")

        pylab.subplot(3,1,2)
        pylab.plot(self.vec)
        pylab.title("proposition")

        pylab.subplot(3,1,3)
        y, x, _ = pylab.hist(self.vec, bins, density=True, lw=0.5, ec="k")
        M1 = max(y)

        # this normalisation is an approximation/hack
        pylab.plot(self.Xtarget, self.Ytarget/ (max(self.Ytarget)/M1), "-ro")
        pylab.title("simulated (blue) and target (red) distributions")

    def check(self, bins=60):
        y, x, _ = pylab.hist(self.vec, bins, density=True, lw=0.5, ec="k")
        M1 = max(y)

        # this normalisation is an approximation/hack
        pylab.plot(self.Xtarget, self.Ytarget/ (max(self.Ytarget)/M1), "-ro")
        pylab.title("simulated (blue) and target (red) distributions")
