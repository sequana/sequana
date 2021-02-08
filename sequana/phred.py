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
"""Manipulate phred quality of reads

FastQ quality are stored as characters. The phred scales indicates the range of characters.

In general, characters goes from ! to ~ that is from 33 to 126 in an ascii
table. This convention starts at 33 because characters before ! may cause trouble 
(e.g. white spaces). This scale is the Sanger scale. There are 2 other scales
that could be used ranging from 59 to 126 (illumina 1) and from 64 to 126 (illumina 1.3+).



So, here are the offset to use:

============== ============ ===============
 Name           offset       Numeric range
============== ============ ===============
 Sanger         33            0 to 93
 Solexa         64            -5 to 62
 illumina1.3+   64            0 to 62
============== ============ ===============

:reference: http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2847217/


Even though dedicated tools would have better performances, we provide a set
of convenient functions. An example is provided here below to plot the quality
corresponding to a character string extracted from a FastQ read.


In this example, we use :class:`Quality` class where the default offset is 33
(Sanger). We compare the quality for another offset 

.. plot::
    :include-source:

    from sequana import phred

    from sequana.phred import Quality
    q = Quality('BCCFFFFFHHHHHIIJJJJJJIIJJJJJJJJFH')
    q.plot()
    q.offset = 64
    q.plot()
    from pylab import legend
    legend(loc="best")


"""
from sequana.lazy import numpy as np
from sequana.lazy import pylab

import colorlog
logger = colorlog.getLogger(__name__)



quality = r"""!"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOP"""
quality += r"""QRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~"""


from math import log10


__all__ = ['Quality', "proba_to_quality_sanger", "quality_to_proba_sanger"]


def proba_to_quality_sanger(pe):
    """A value between 0 and 93

    :param pe: the probability of error.
    :return: Q is the quality score.

    - a high probability of error (0.99) gives Q=0
    - q low proba of errors (0.05) gives Q = 13
    - q low proba of errors (0.01) gives Q = 20

    """
    if pe > 1:
        pe = 1
    if pe < 1e-90:
        pe = 1e-90
    Qs = -10 * log10(pe)
    if Qs > 93:
        Qs = 93
    return Qs

def quality_to_proba_sanger(quality):
    """Quality to probability (Sanger)"""
    return 10**(quality/-10.)


def proba_to_quality_solexa(pe):
    """prior v1.3 (ref: wikipedia
    https://en.wikipedia.org/wiki/FASTQ_format
    """
    if pe > 1:
        pe = 1
        return -5

    if pe <1e-90:
        pe = 1e-90
    Qs = -10 * log10(pe/(1-pe))
    if Qs > 62:
        Qs = 62
    if Qs < -5:
        Qs = -5
    return Qs


def quality_solexa_to_quality_sanger(qual):
    return 10 * log10(10**(qual/10.) + 1 )


def quality_sanger_to_quality_solexa(qual):
    """


    """
    return 10 * log10(10**(qual/10.) - 1 )


def ascii_to_quality(character, phred=33):
    """ASCII to Quality conversion

    :param int phred: offset (defaults to 33)

    ::

        >>> ascii_to_quality("!")
        0
        >>> ascii_to_quality("~")
        93
    """
    return ord(character) - phred


def quality_to_ascii(quality, phred=33):
    """Quality to ASCII conversion 

    :param int phred: offset (defaults to 33)

    ::

         >>> quality_to_ascii(65)
         b

    """
    return chr(quality + phred)


class Quality(object):
    """Phred quality

    .. plot::
        :include-source:

        >>> from sequana.phred import Quality
        >>> q = Quality('BCCFFFFFHHHHHIIJJJJJJIIJJJJJJJJFH')
        >>> q.plot()

    You can access to the quality as a list using the :attr:`quality` attribute
    and the mean quality from the :attr:`mean_quality` attribute.

    """
    def __init__(self, seq, offset=33):
        self.seq = seq
        self.offset = offset

    def _get_quality(self):
        # bytearray conversion is required since ord() function
        # does not handle bytes from py3
        return [x - self.offset for x in bytearray(self.seq, 'utf-8')]
    quality = property(_get_quality, doc="phred string into quality list")

    def _get_mean_quality(self):
        return np.mean(self.quality)
    mean_quality = property(_get_mean_quality, doc="return mean quality")

    def plot(self, fontsize=16):
        """plot quality versus base position"""
        pylab.plot(self.quality, label="offset: %s" % self.offset)
        pylab.xlabel('base position', fontsize=fontsize)
        pylab.ylabel('Quality per base', fontsize=fontsize)
        pylab.grid(True)
        # ylim set autoscale to off so if we want to call this function  several
        # times, we must reset autoscale to on before calling ylim
        pylab.autoscale()
        limits = pylab.ylim()
        pylab.ylim(max(0,limits[0]-1), limits[1]+1)



# this should qlso be correct for Illumina 1.8+
class QualitySanger(Quality):
    """Specialised :class:`Quality` class for Sanger case"""
    def __init__(self, seq):
        super(QualitySanger, self).__init__(seq, offset=33)


class QualitySolexa(Quality):
    """Specialised :class:`Quality` class for Solexa case"""
    def __init__(self, seq):
        super(QualitySolexa, self).__init__(seq, offset=64)


















