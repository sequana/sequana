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
from sequana.lazy import pandas as pd
from sequana.lazy import pylab
from pylab import plot

import colorlog
logger = colorlog.getLogger(__name__)




class CNVnator(object):
    """Reader of the CNVnator output file.



    """
    def __init__(self, filename):

        self.filename = filename


        self.df = pd.read_csv(filename, sep="\t", header=None)
        self.df.columns = ["type", "label", "size", "3","4","5","6","7","8"]
        name = self.df['label'].apply(lambda x: x.rsplit(":",1)[0])
        start = self.df['label'].apply(lambda x: x.rsplit(":",1)[1].split("-")[0])
        end = self.df['label'].apply(lambda x: x.rsplit(":",1)[1].split("-")[1])
        self.df['start'] = start.astype(int)
        self.df['end'] = end.astype(int)
        self.df['name'] = name

    def plot(self, chr_name, x1=None, x2=None, Y=20):

        df = self.df.query("name == @chr_name")
        for _, item in df.iterrows():
            if item['type'] == "deletion":
                plot([item.start, item.end], [-1,-1], "r-", label="deletion")
            else:
                plot([item.start, item.end], [Y, Y], "b-", label="duplication")
        pylab.legend()

    def hist_event_size(self, bins=20):
        self.df['size'].hist(ec="k", bins=bins, lw=1)

