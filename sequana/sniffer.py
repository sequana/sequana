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
"""Sniffer"""
from sequana.bamtools import is_bam, is_sam, is_cram
from sequana.fastq import is_fastq
from sequana.fasta import is_fasta

import colorlog
logger = colorlog.getLogger(__name__)



def sniffer(filename):

    try:
        if is_sam(filename): return "SAM"
    except:
        pass

    try:
        if is_bam(filename): return "BAM"
    except:
        pass

    try:
        if is_cram(filename): return "CRAM"
    except:
        pass

    try:
        if is_fastq(filename): return "FASTQ"
    except:
        pass

    try:
        if is_fasta(filename): return "FASTA"
    except:
        pass

