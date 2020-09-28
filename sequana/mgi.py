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


from sequana.lazy import pandas as pd
from sequana.lazy import pylab
from sequana.lazy import numpy as np


class MGI():

    def __init__(self, filename):
        self.filename = filename
        self.df = pd.read_csv(self.filename, sep='\t', comment="#", header=None)

        self.metadata = {}
        with open(self.filename, "r") as fin:
            line = fin.readline()
            while line.startswith('#'):
                key, value = line.split(maxsplit=1)
                key = key[1:]  # skip the # character
                try:
                    self.metadata[key] = int(value.strip())
                except:
                    self.metadata[key] = value.strip()
                if key.endswith("%"):
                    try:
                        self.metadata[key] = float(value.strip())
                    except: #pragma: no cover
                        pass
                line = fin.readline()
        try:
            ncounts = self.metadata['N_Count']
            self.metadata['N_Count'] =  int(ncounts.split('\t')[0])
            self.metadata['N_Count%'] =  float(ncounts.split('\t')[1]  )
        except Exception as err: #pragma no cover
            print(err)
            pass
        self.df.columns = ['Pos'] + self.metadata['Pos'].split('\t')
        del self.metadata['Pos']
        self.df.set_index('Pos', inplace=True)


    def plot_acgt(self):
        ACGTN = self.df[['A', 'C', 'G', 'T', 'N']]
        N = ACGTN.sum(axis=1).iloc[0]
        assert N == self.metadata['ReadNum']
        for this in 'ACGTN':
            pylab.plot(self.df[this] / N * 100, label=this)


    def boxplot_quality(self, color_line='r', bgcolor='grey', color='yellow', lw=4, 
            hold=False, ax=None):


        quality = self.df[[str(x) for x in range(42)]]  # not sure why we have phred score from 0 to 41
        N = self.metadata['ReadNum']
        proba = quality / N

        self.xmax = 150
        xmax = self.xmax + 1
        if ax:
            pylab.sca(ax) # pragma no cover
        pylab.fill_between([0,xmax], [0,0], [20,20], color='red', alpha=0.3)
        pylab.fill_between([0,xmax], [20,20], [30,30], color='orange', alpha=0.3)
        pylab.fill_between([0,xmax], [30,30], [41,41], color='green', alpha=0.3)


        X = []
        Q = []
        S = []
        for pos in range(1, 151):
            qualities = [((int(k)+1)*v) for k,v in quality.loc[pos].items()]
            mean_quality = sum(qualities) / N
            X.append(pos)
            Q.append(mean_quality)
            proba = quality.loc[pos] / N

            std = pylab.sqrt(sum([(x-mean_quality)**2 * y for x, y in zip(range(42), proba)]))
            S.append(std)

        print(len(X))
        print(len(Q))
        print(len(S))

        Q = np.array(Q)
        X = np.array(X)
        S = np.array(S)
        pylab.fill_between(X, Q+S, Q-S, 
            color=color, interpolate=False)

        pylab.plot(X, Q, color=color_line, lw=lw)
        pylab.ylim([0, 41])
        pylab.xlim([0, self.xmax+1])
        pylab.title("Quality scores across all bases")
        pylab.xlabel("Position in read (bp)")
        pylab.ylabel("Quality")
        pylab.grid(axis='x')















