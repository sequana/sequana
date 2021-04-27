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



class Boxplot(object):
    """Used to plot boxplot of fastq quality a la fastqc"""
    
    def __init__(self, data):

        # if data is a dataframe, keep it else, transform to dataframe
        try:
            self.df = pd.DataFrame(data)
        except:
            self.df = data

        self.xmax = self.df.shape[1]
        self.X =  None

    def plot(self, color_line='r', bgcolor='grey', color='yellow', lw=4, 
            hold=False, ax=None):

        xmax = self.xmax + 1
        if ax:
            pylab.sca(ax)

        if self.X is None:
            X = range(1, self.xmax + 1)

        pylab.fill_between(X, 
            self.df.mean()+self.df.std(), 
            self.df.mean()-self.df.std(), 
            color=color, interpolate=False)

        pylab.plot(X, self.df.mean(), color=color_line, lw=lw)

        if self.df.mean().mean()<40: #illumina
            pylab.fill_between([0,xmax], [0,0], [20,20], color='red', alpha=0.3)
            pylab.fill_between([0,xmax], [20,20], [30,30], color='orange', alpha=0.3)
            pylab.fill_between([0,xmax], [30,30], [41,41], color='green', alpha=0.3)
            pylab.ylim([0, 41])
        else:
            pylab.fill_between([0,xmax], [0,0], [20,20], color='red', alpha=0.3)
            pylab.fill_between([0,xmax], [20,20], [30,30], color='orange', alpha=0.3)
            pylab.fill_between([0,xmax], [30,30], [101,101], color='green', alpha=0.3)
            pylab.ylim([0, 101])

        pylab.xlim([0, self.xmax+1])
        pylab.title("Quality scores across all bases")
        pylab.xlabel("Position in read (bp)")
        pylab.ylabel("Quality")
        pylab.grid(axis='x')
