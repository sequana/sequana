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
import itertools

import colorlog
logger = colorlog.getLogger(__name__)



def build_kmer(length=6, letters='CG'):
    """Return list of kmer of given length based on a set of letters

    :return: list of kmers
    """
    # build permutations of CG letters with a sequence of given lengths
    # TODO include N other letters
    combos = list(itertools.product(letters, repeat=length))
    return ["".join(this) for this in combos]


def get_kmer(sequence, k=7):
    """Given a sequence, return consecutive kmers


    :return: iterator of kmers

    """
    for i in range(0, len(sequence)-k+1):
        yield sequence[i:i+k]
