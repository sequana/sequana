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
import math

from sequana.lazy import numpy as np
from sequana.lazy import pandas as pd

import colorlog

logger = colorlog.getLogger(__name__)


__all__ = ["Coverage"]


class Coverage(object):
    r"""Utilities related to Lander and Waterman theory

    We denote :math:`G` the genome length in nucleotides and :math:`L` the read
    length in nucleotides. These two numbers are in principle well defined since
    :math:`G` is defined by biology and :math:`L` by the sequencing machine.

    The total number of reads sequenced during an experiment is denoted
    :math:`N`. Therefore the total number of nucleotides is simply :math:`NL`.

    The depth of coverage (DOC) at a given nucleotide position is the number of
    times that a nucleotide is covered by a mapped read.

    The theoretical fold-coverage is defined as :

    .. math:: a = NL/G

    that is the average number of times each nucleotide is expected to be sequenced (in the
    whole genome). The fold-coverage is often denoted :math:`aX` (e.g., 50X).

    In the :class:`Coverage` class, :math:`G` and :math:`N` are fixed at
    the beginning. Then, if one changes :math:`a`, then :math:`N` is updated and
    vice-versa so that the relation :math:`a=NL/G` is always true::

        >>> cover = Coverage(G=1000000, L=100)
        >>> cover.N = 100000    # number of reads
        >>> cover.a             # What is the mean coverage
        10
        >>> cover.a = 50
        >>> cover.N
        500000

    From the equation aforementionned, and assuming the reads are uniformly
    distributed, we can answer a few interesting questions using probabilities.

    In each chromosome, a read of length :math:`L` could start at any position
    (except the last position L-1). So in a genome :math:`G` with :math:`n_c`
    chromosomes, there are :math:`G - n_c (L-1)` possible starting positions.
    In general  :math:`G >> n_c (L-1)` so the probability that one of the
    :math:`N` read starts at any specific nucleotide is :math:`N/G`.

    The probability that a read (of length :math:`L`) covers a given
    position is :math:`L/G`. The probability of **not** covering that location
    is :math:`1-L/G`. For :math:`N` fragments, we obtain the probability
    :math:`\left(1-L/G\right)^N`. So, the probability of covering a given
    location with at least one read is :

    .. math:: P = 1 - \left(1- \frac{L}{G}\right)^N

    Since in general, N>>1, we have:

    .. math::  P = 1- \exp^{-NL/G}

    From this equation, we can derive the fold-coverage required to have
    e.g., :math:`E=99\%` of the genome covered:

    .. math:: a = log(-1/(E-1)

    .. indeed e^x = lim (1+x/N)^N when N goes to infinite so for x -> -x.N
        e^-Nx = lim (1-x)^N when N goes to infinite

    equivalent to

    .. math:: a = -log(1-E)

    The method :meth:`get_required_coverage` uses this equation. However, for
    numerical reason, one should not provide :math:`E` as an argument but (1-E).
    See :meth:`get_required_coverage`

    Other information can also be derived using the methods
    :meth:`get_mean_number_contigs`, :meth:`get_mean_contig_length`,
    :meth:`get_mean_contig_length`.

    .. seealso:: :meth:`get_table` that provides a summary of all these
        quantities for a range of coverage.

    :reference: http://www.math.ucsd.edu/~gptesler/283/FA14/slides/shotgun_14-handout.pdf

    """

    def __init__(self, N=None, L=None, G=None, a=None):
        self._a = None  # coverage
        self._L = L  # length of each read
        self._N = N  # number of reads
        self._G = G  # length target genome

    def __repr__(self):
        if self._N is None:
            N = "undefined"
        else:
            N = self._N
        if self._a is None:
            a = "undefined"
        else:
            a = self._a
        return "Coverage(N=%s, L=%s, G=%s, a=%s) " % (N, self.L, self.G, a)

    def _get_a(self):
        return self._N * self._L / float(self._G)

    def _set_a(self, value):
        assert value > 0
        self._a = value
        self._N = self._a * self._G / float(self._L)

    a = property(_get_a, _set_a, doc="coverage defined as NL/G")

    def _get_L(self):
        return self._L

    def _set_L(self, value):
        assert value > 0
        self._L = value

    L = property(_get_L, _set_L, doc="length of the reads")

    def _get_N(self):
        return self._N

    def _set_N(self, value):
        assert value > 0
        self._N = value

    N = property(_get_N, _set_N, doc="number of reads defined as aG/L")

    def _get_G(self):
        return self._G

    def _set_G(self, value):
        assert value > 0
        self._G = value

    G = property(_get_G, _set_G, doc="genome length")

    def get_required_coverage(self, M=0.01):
        """Return the required coverage to ensure the genome is covered

        A general question is what should be the coverage to make sure
        that e.g. E=99% of the genome is covered by at least a read.

        The answer is:

        .. math:: \log^{-1/(E-1)}

        This equation is correct but have a limitation due to floating precision.
        If one provides E=0.99, the answer is 4.6 but we are limited to a
        maximum coverage of about 36 when one provides E=0.9999999999999999
        after which E is rounded to 1 on most computers. Besides, it is no
        convenient to enter all those numbers. A scientific notation would be better but
        requires to work with :math:`M=1-E` instead of :math:`E`.

        .. math:: \log^{-1/ - M}

        So instead of asking the question what is the
        requested fold coverage to have 99% of the genome covered, we ask the question what
        is the requested fold coverage to have 1% of the genome not covered.
        This allows us to use :math:`M` values as low as 1e-300 that is a fold coverage
        as high as 690.


        :param float M: this is the fraction of the genome not covered by
            any reads (e.g. 0.01 for 1%). See note above.
        :return: the required fold coverage

        .. plot::

            import pylab
            from sequana import Coverage
            cover = Coverage()
            misses = np.array([1e-1, 1e-2, 1e-3, 1e-4,1e-5,1e-6])
            required_coverage = cover.get_required_coverage(misses)
            pylab.semilogx(misses, required_coverage, 'o-')
            pylab.ylabel("Required coverage", fontsize=16)
            pylab.xlabel("Uncovered genome", fontsize=16)
            pylab.grid()

        # The inverse equation is required fold coverage = [log(-1/(E - 1))]
        """
        # What should be the fold coverage to have 99% of the genome sequenced ?
        # It is the same question as equating 1-e^{-(NL/G}) == 0.99, we need NL/G = 4.6
        if isinstance(M, float) or isinstance(M, int):
            assert M < 1
            assert M >= 0
        else:
            M = np.array(M)
        # Here we do not use log(-1/(E-1)) but log(-1/(1-E-1)) to allow
        # for using float down to 1e-300 since 0.999999999999999 == 1
        return np.log(-1 / (-M))

    def get_mean_number_contigs(self):
        r"""Expected number of contigs

        A binomial distribution with parameters :math:`N` and :math:`p`

        .. math:: \left(\frac{aG}{L}\right) \exp^{-a}

        """
        return self.G / float(self.L) * self.a * math.exp(-self.a)

    def get_mean_contig_length(self):
        """Expected length of the contigs

        .. math:: \\frac{e^a-1)L}{a}

        """
        return (math.exp(self.a) - 1) * self.L / self.a

    def get_mean_reads_per_contig(self):
        r"""Expected number of reads per contig

        Number of reads divided by expected number of contigs:

        .. math:: N / (N\exp^{-a}) = e^a


        """
        return math.exp(self.a)

    def get_percent_genome_sequenced(self):
        """Return percent of the genome covered

        .. math:: 100 (1-\exp{-a})
        """
        return 100 * (1 - math.exp(-self.a))

    def get_summary(self):
        """Return a summary (dictionary) for the current fold coverage :attr:`a`"""
        data = {}
        data["coverage"] = self.a
        data["reads"] = self.N
        data["nucleotides"] = self.a * self.G
        data["% genome sequenced"] = self.get_percent_genome_sequenced()
        data["mean number of contigs"] = self.get_mean_number_contigs()
        data["mean contig length"] = self.get_mean_contig_length()
        data["mean reads per contig"] = self.get_mean_reads_per_contig()
        return data

    def get_table(self, coverage=None):
        """Return a summary dataframe for a set of fold coverage

        :param list coverage: if None, coverage list starts at 0.5 and ends at
            10 with a step of 0.5
        """

        if coverage is None:
            X = np.arange(0.5, 10.5, 0.5)
        else:
            X = coverage

        a_buf = self.a

        results = []
        for this in X:
            self.a = this
            results.append(self.get_summary())

        df = pd.DataFrame.from_records(results)
        df.set_index("coverage", inplace=True)
        return df
