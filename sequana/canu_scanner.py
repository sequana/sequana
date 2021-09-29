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

import os
import glob
from sequana import version
from sequana.lazy import pandas as pd
from sequana.lazy import pylab

import colorlog
logger = colorlog.getLogger(__name__)




class CanuScanner():
    """


    """
    def __init__(self, path="."):
        """

        """
        self.path = path

        self.data = {}
        self.data["version"] = version
        self.data["tool"] = "sequana"
        self.data["module"] = "canu_scanner"
        self.data["correction"] = {}
        self.data["trimming"] = {}
        self.data["assembly"] = {}

    def getfile(self, filename):
        filenames = glob.glob(self.path + os.sep + filename)
        assert len(filenames) == 1
        return filenames[0]

    ######################################################################""# CORRECTION 

    def scan_correction(self):
        """

        """
        filename = self.getfile("correction/*.gkpStore/load.dat")
        with open(filename, "r") as fin:
            data = fin.read().split("\n")

        # Get kept reads
        bp = int(data[2].split()[2])
        reads = int(data[2].split()[1])
        self.data['correction']["readsLoaded"] = {"reads":reads, "bp":bp}

        # Get skipped reads
        bp = int(data[2].split()[4])
        reads = int(data[2].split()[3])
        self.data['correction']["readsSkipped"] = {"reads":reads, "bp":bp}

    def hist_read_length(self, bins=100, fontsize=16):
        """

        """
        filename = self.getfile("correction/*gkpStore/reads.txt")

        df = pd.read_csv(filename, header=None, sep="\t")
        df[2].hist(bins=bins)
        df.columns = ["ID", 1, "read_length", 3, 4]
        pylab.xlabel("Read length", fontsize=fontsize)
        pylab.ylabel("Number of reads", fontsize=fontsize)
        pylab.xlim([0, pylab.xlim()[1]])
        return df

    def plot_kmer(self, bins=100):
        filename = self.getfile("correction/0-mercounts/*.ms16.histogram")

        df = pd.read_csv(filename, header=None, sep="\t")
        df.columns = ["kmer", "count", "X", "Y"]

        # Save some data
        self.data['correction']['largest mercount'] = list(df['kmer'])[-1]
        self.data['correction']['unique mers'] = df['count'][0]
        self.data['correction']['distinc mers'] = df['count'].sum()
        self.data['correction'][""] = sum(df.kmer * df['count'])

        # X is just df['count'].cumsum() / df['count'].sum() (distinct kmer)
        # Y is (df['kmer']*df['count']).cumsum() / (df['kmer']*df['count()).sum()
        # that is total kmer 
        pylab.plot(df.X, df.Y, label="distinct vs total")
        pylab.grid()
        pylab.legend()
        return df

    def set_overlap_filtering(self):
        filename = self.getfile("correction/2-correction/*.globalScores.stats")

        with open(filename, "r") as fin:
            data = fin.read()
        self.data['correction']["overlap filtering"] = data

    def set_read_correction(self):
        filename = self.getfile("correction/2-correction/*.correction.summary")

        with open(filename, "r") as fin:
            data = fin.read()
        self.data['correction']["read correction"] = data

    def plot_correction_check1(self, alpha=0.5):
        try:
            tn = pd.read_csv(
                self.getfile("correction/2-correction/*.estimate.tn.log"),
                sep="\t", header=None, usecols=[0, 1, 2, 3])
            pylab.plot(tn[1], tn[3], "x", color="purple", label="TN", alpha=alpha)
        except:
            pass

        try:
            fn = pd.read_csv(
                self.getfile("correction/2-correction/*.estimate.fn.log"),
                sep="\t", header=None, usecols=[0, 1, 2, 3])
            pylab.plot(fn[1], fn[3], "x", color="green", label="FN", alpha=alpha)
        except:
            pass
    
        try:
            fp = pd.read_csv(
                self.getfile("correction/2-correction/*.estimate.fp.log"),
                sep="\t", header=None, usecols=[0, 1, 2, 3])
            pylab.plot(fp[1], fp[3], "x", color="cyan", label="FP", alpha=alpha)
        except:
            pass

        try:
            tp = pd.read_csv(
                self.getfile("correction/2-correction/*.estimate.tp.log"),
                sep="\t", header=None, usecols=[0, 1, 2, 3])
            pylab.plot(tp[1], tp[3], "x", color="orange", label="TP", alpha=alpha)
        except:
            pass

        pylab.xlabel("original read length")
        pylab.ylabel("corrected read length")
        pylab.legend()
        pylab.grid()
        pylab.xlim(0, pylab.xlim()[1])
        pylab.ylim(0, pylab.xlim()[1])


        caption = """ Scatter plot of the original read length (X axis) against
the expected corrected read length (Y axis). Colors show a comparison of the
simple filter (which doesn't use overlaps) to the expensive filter (which does).
A large green triangle (false negatives) hints that there could be abnormally
low quality regions in the reads. """ # from canu report

    def hist_read_length2(self, fontsize=16):
        df = pd.read_csv(
            self.getfile("correction/2-correction/*.original-expected-corrected-length.dat"),
            sep="\t", header=None)
        pylab.clf()
        df[1].hist(bins=100, alpha=0.5, density=True, label="original")
        df[2].hist(bins=100, alpha=0.5, density=True, label="expected")
        df[3].hist(bins=100, alpha=0.5, density=True, label="corrected")
        pylab.legend()
        pylab.xlabel("read length", fontsize=fontsize)
        pylab.ylabel("number of reads ", fontsize=fontsize)
        return df

    ######################################################################""# CORRECTION 

    def scan_trimming(self):
        filename = self.getfile("trimming/*.gkpStore/load.dat")
        with open(filename, "r") as fin:
            data = fin.read().split("\n")

        # Get kept reads
        bp = int(data[2].split()[2])
        reads = int(data[2].split()[1])
        self.data['trimming']["readsLoaded"] = {"reads":reads, "bp":bp}

        # Get skipped reads
        bp = int(data[2].split()[4])
        reads = int(data[2].split()[3])
        self.data['trimming']["readsSkipped"] = {"reads":reads, "bp":bp}

    def hist_trimming_read_length(self, bins=100, fontsize=16):
        """

        """
        filename = self.getfile("trimming/*gkpStore/reads.txt")

        df = pd.read_csv(filename, header=None, sep="\t")
        df[2].hist(bins=bins)
        df.columns = ["ID", 1, "read_length", 3, 4]
        pylab.xlabel("Read length", fontsize=fontsize)
        pylab.ylabel("Number of reads", fontsize=fontsize)
        pylab.xlim([0, pylab.xlim()[1]])
        return df


