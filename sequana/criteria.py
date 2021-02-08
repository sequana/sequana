# -*- coding: utf-8 -*-
#
#  This file is part of Sequana software
#
#  Copyright (c) 2016 - Sequana Development Team
#
#  File author(s):
#      Thomas Cokelaer <thomas.cokelaer@pasteur.fr>
#      Dimitri Desvillechabrol <dimitri.desvillechabrol@pasteur.fr>, 
#          <d.desvillechabrol@gmail.com>
#
#  Distributed under the terms of the 3-clause BSD license.
#  The full license is in the LICENSE file, distributed with this software.
#
#  website: https://github.com/sequana/sequana
#  documentation: http://sequana.readthedocs.io
#
##############################################################################
"""Akaike and other criteria

"""
import math

import colorlog
logger = colorlog.getLogger(__name__)




__all__ = ['AIC', 'AICc', 'BIC']


def AIC(L, k, logL=False):
    r"""Return Akaike information criterion (AIC)

    :param float L: maximised value of the likelihood function
    :param int k: number of parameters
    :param bool logL: L is the log likelihood.

    Suppose that we have a statistical model of some data, from which we computed
    its likelihood function and let :math:`k` be the number of parameters in the model
    (i.e. degrees of freedom). Then the AIC value is:

    :math:`\mathrm{AIC} = 2k - 2\ln(L)`

    Given a set of candidate models for the data, the preferred model is the one
    with the minimum AIC value. Hence AIC rewards goodness of fit (as assessed
    by the likelihood function), but it also includes a penalty that is an
    increasing function of the number of estimated parameters. The penalty
    discourages overfitting.

    Suppose that there are R candidate models AIC1, AIC2, AIC3, AICR.
    Let AICmin be the minimum of those values. Then, exp((AICmin - AICi)/2)
    can be interpreted as the relative probability that the ith model
    minimizes the (estimated) information loss.

    Suppose that there are three candidate models, whose AIC values are 100,
    102, and 110. Then the second model is exp((100 - 102)/2) = 0.368 times
    as probable as the first model to minimize the information loss. Similarly,
    the third model is exp((100 - 110)/2) = 0.007 times as probable as
    the first model, which can therefore be discarded.

    With the remaining two models, we can (1) gather more data, (2) conclude
    that the data is insufficient to support selecting one model from among
    the first two (3) take a weighted average of the first two models,
    with weights 1 and 0.368.

    The quantity `exp((AIC_{min} - AIC_i)/2)` is the relative likelihood of model i.

    If all the models in the candidate set have the same number of parameters,
    then using AIC might at first appear to be very similar to using the
    likelihood-ratio test. There are, however, important distinctions.
    In particular, the likelihood-ratio test is valid only for nested models,
    whereas AIC (and AICc) has no such restriction.

    :Reference: Burnham, K. P.; Anderson, D. R. (2002), Model Selection and
        Multimodel Inference: A Practical Information-Theoretic Approach (2nd ed.),
        Springer-Verlag, ISBN 0-387-95364-7.
    """
    if logL is True:
        return 2 * k + 2 * L
    else:
        return 2 * k -2 * math.log(L)


def AICc(L, k, n, logL=False):
    r"""AICc criteria

    :param float L: maximised value of the likelihood function
    :param int k: number of parameters
    :param int n: sample size
    :param bool logL: L is the log likelihood.

    AIC with a correction for finite sample sizes.
    The formula for AICc depends upon the statistical model.
    Assuming that the model is univariate, linear, and has normally-distributed
    residuals (conditional upon regressors), the formula for AICc is as follows:

    AICc is essentially AIC with a greater penalty for extra parameters.
    Using AIC, instead of AICc, when n is not many times larger than k2, increases
    the probability of selecting models that have too many parameters, i.e. of
    overfitting. The probability of AIC overfitting can be substantial, in some cases.

    """
    res = AIC(L, k, logL=logL) + 2*k*(k+1.) / (n-k-1.)
    return res



def BIC(L, k, n, logL=False):
    r"""Bayesian information criterion

    :param float L: maximised value of the likelihood function
    :param int k: number of parameters
    :param int n: sample size
    :param bool logL: L is the log likelihood.

    Given any two estimated models, the model with the lower value of BIC is the one to be preferred.
    """

    if logL is False:
        res = -2 * math.log(L) + k * (math.log(n) - math.log(2 * math.pi))
        # For large n
        #res = -2 * math.log(L) + k * math.log(n)
        return res
    else:
        res = 2 * logL + k * (math.log(n) - math.log(2 * math.pi))
        return res







