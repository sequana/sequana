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

from sequana.lazy import pylab
from sequana.utils.pandas import PandasReader

import colorlog

from collections import Counter

logger = colorlog.getLogger(__name__)


__all__ = ["Homer"]


class Homer:
    def __init__(self, filename):

        header = open(filename).readline().strip().split("\t")[1:]

        self.df = PandasReader(filename, sep="\t", skiprows=1, header=None).df
        self.df.columns = ["ID"] + header
        self.df.fillna("NA", inplace=True)

    def pie_annotation(self, wedgeprops={"ec": "k"}, **kwargs):

        counts = Counter(self.df.Annotation)
        labels = [f"{k.split()[0]} ({str(v)})" for k, v in counts.items()]
        pylab.pie(counts.values(), labels=labels, wedgeprops=wedgeprops, **kwargs)
