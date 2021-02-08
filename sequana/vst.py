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
"""Simple VST transformation"""
from sequana.lazy import numpy as np
import colorlog
logger = colorlog.getLogger(__name__)


__all__ = ['VST']

class VST():
    """

    ::
    
        v = VST()
        v.anscombe(X)

    Reference: M. Makitalo and A. Foi, "Optimal inversion of the
        generalized Anscombe transformation for Poisson-Gaussian noise",
        IEEE Trans. Image Process., doi:10.1109/TIP.2012.2202675

    """
    def __init__(self):
        pass

    @staticmethod
    def anscombe(x):
        r"""Compute the anscombe variance stabilizing transform.

        :param x: noisy Poisson-distributed data
        :return: data with variance approximately equal to 1.

        Reference: Anscombe, F. J. (1948), "The transformation of Poisson,
            binomial and negative-binomial data", Biometrika 35 (3-4): 246-254

        For Poisson distribution, the mean and variance are not independent. The
        anscombe transform aims at transforming the data so that the variance
        is about 1 for large enough mean; For mean zero, the varaince is still
        zero. So, it transform Poisson data to approximately Gaussian data with
        mean :math:`\sqrt{x+3/8} - 1/(4m^{1/2})`
        """

        #if np.mean(x) <4:
        #    logger.warning("Mean of input data below 4")
        try:
            # If a dataframe, we do not want to change it
            return 2.0 * np.sqrt(x + 3.0/8.0)
        except:
            return 2.0 * np.sqrt(np.array(x) + 3.0/8.0)

    @staticmethod
    def inverse_anscombe(x):
        """Compute the inverse transform

        This uses an approximation of the exact  unbiased inverse.
        Reference: Makitalo, M., & Foi, A. (2011). A closed-form
        approximation of the exact unbiased inverse of the Anscombe
        variance-stabilizing transformation. Image Processing.
        """
        # should be positive
        # should wotk for dataframe and numpy arrays 
        return (1.0/4.0 * np.power(x, 2) +
                1.0/4.0 * np.sqrt(3.0/2.0) * np.power(x, -1.0) -
                11.0/8.0 * np.power(x, -2.0) + 
                5.0/8.0 * np.sqrt(3.0/2.0) * np.power(x, -3.0) - 1.0 / 8.0)

    @staticmethod
    def generalized_anscombe(x, mu, sigma, gain=1.0):
        """Compute the generalized anscombe variance stabilizing transform

        Data should be a mixture of poisson and gaussian noise.

        The input signal  z  is assumed to follow the Poisson-Gaussian noise
        model::
    
            x = gain * p + n

        where gain is the camera gain and mu and sigma are the read noise
        mean and standard deviation. X should contain only positive values.  
        Negative values are ignored. Biased for low counts
        """
        try:
            # If a dataframe, we do not want to change it
            y = gain*x + (gain**2) * 3.0/8.0 + sigma**2 - gain*mu
            return (2.0 / gain) * np.sqrt(np.maximum(y, 0.0))
        except:
            x = np.array(x)
            y = gain*x + (gain**2) * 3.0/8.0 + sigma**2 - gain*mu
            return (2.0 / gain) * np.sqrt(np.maximum(y, 0.0))

