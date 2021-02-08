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
"""
Python script to filter a VCF file
"""
import sys

from sequana.lazy import vcf

import colorlog
logger = colorlog.getLogger(__name__)



__all__ = ["VCFBase"]


class VCFBase(vcf.Reader):
    """Base class for VCF files

    Read an existing file as follows::

        from sequana.vcf_filter import VCFBase
        v = VCFBase("filename.vcf")

    You can get the number of variants::

        len(v)

    the version and source of the VCF creator (if provided)::

        v.version, v.source

    and you can easily iterate through the variants::

        for variant in vcf:
            print(variant)

    note that if you iterate again, you will get nothing. You will need to
    rewind the cursor::

        vcf.rewind()

    You also get lots of extra information inherited from the vcf.Reader 

    """
    def __init__(self, filename, verbose=True, **kwargs):
        """.. rubric:: constructor

        :param str filename: a vcf file.
        :param kwargs: any arguments accepted by vcf.Reader

        """
        try:
            self.filename = filename
            filin = open(filename, "r")
            vcf.Reader.__init__(self, fsock=filin, **kwargs)
            self._get_start_index()
        except FileNotFoundError as e:
            logger.error(
                "FileNotFoundError({0}): {1}".format(e.errno, e.strerror)
            )
            raise FileNotFoundError

        if verbose:
            print("Found VCF version {}".format(self.version))

    def rewind(self):
        """ Rewind the reader
        """
        self._reader.seek(self._start_index)
        self.reader = (line.strip() for line in self._reader if line.strip())

    def __len__(self):
        self.rewind()
        i = 0
        for line in self:
            i += 1
        self.rewind()
        return i

    def _get_start_index(self):
        self._reader.seek(0)
        for line in iter(self._reader.readline, ''):
            if line.startswith("#"):
                self._start_index = self._reader.tell()
            else:
                self.rewind()
                break

    def _get_version(self):
        fileformat = self.metadata['fileformat']
        if fileformat == 'VCFv4.1':
            return "4.1"
        elif fileformat == "VCFv4.2":
            return "4.2"
        else:
            return fileformat
    version = property(_get_version)

    def _get_source(self):
        if "source" in self.metadata:
            return self.metadata['source'][0]
        else:
            return "undefined"
    source = property(_get_source)


