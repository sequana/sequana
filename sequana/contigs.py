#
#  This file is part of Sequana software
#
#  Copyright (c) 2016-2022 - Sequana Development Team
#
#  Distributed under the terms of the 3-clause BSD license.
#  The full license is in the LICENSE file, distributed with this software.
#
#  website: https://github.com/sequana/sequana
#  documentation: http://sequana.readthedocs.io
#
##############################################################################
from sequana.lazy import pylab
from sequana.lazy import pandas as pd
from sequana.lazy import numpy as np

from sequana import tools
from sequana import FastA

import colorlog

logger = colorlog.getLogger(__name__)


__all__ = ["ContigsBase", "Contigs"]


class ContigsBase(object):
    """Parent class for contigs data"""

    def __init__(self, filename):
        """.. rubric:: Constructor

        :param filename: input file name
        """
        self.filename = filename
        self.fasta = FastA(filename)

    def get_gc(self, window=100):
        """Return GC content for each contig"""
        data = tools._base_content(self.filename, window, "GC")
        names = self.fasta.names
        lengths = self.fasta.lengths
        GC = [100 * np.nanmean(data[name]) for name in names]
        return GC

    def plot_contig_length_vs_GC(self, alpha=0.5):
        """Plot contig GC content versus contig length

        .. plot::

            from sequana.contigs import Contigs
            from sequana import sequana_data
            filename = sequana_data("test_contigs_spades.fasta")
            ctg = Contigs(filename)
            ctg.plot_contig_length_vs_GC()

        """
        pylab.plot(self.df["length"], self.df["GC"], "o", alpha=alpha)
        pylab.xlabel("contig length (bp)")
        pylab.ylabel("GC (%)")
        pylab.grid(True)
        pylab.ylim([0, 100])
        pylab.xlim(0, self.df["length"].max() + 10)

    def scatter_length_cov_gc(self, min_length=200, min_cov=10, grid=True, logy=False, logx=True):

        """Plot scatter length versus GC content

        :param min_length: add vertical line to indicate possible
            contig length cutoff
        :param min_cov: add horizontal line to indicate possible
            coverage contig cutff
        :param grid: add grid to the plot
        :param logy: set y-axis log scale
        :param logx: set x-axis log scale

        .. plot::

            from sequana import Contigs, sequana_data
            filename = sequana_data("test_contigs_spades.fasta")
            ctg = Contigs(filename)
            ctg.scatter_length_cov_gc()
        """
        if "cov" not in self.df.columns:
            logger.warning("scatter_length_cov_gc required 'cov' coverage column information")
            return
        pylab.clf()
        pylab.scatter(self.df.length, self.df["cov"], c=self.df.GC)
        if logx:
            pylab.semilogx()
        if logy:
            pylab.semilogy()
        pylab.axvline(min_length, lw=2, c="r", ls="--")
        pylab.axhline(min_cov, lw=2, c="r", ls="--")
        pylab.xlabel("contig length")
        pylab.ylabel("contig coverage")
        pylab.colorbar(label="GC")
        if grid:
            pylab.grid(True)


class Contigs(ContigsBase):
    """Utilities for summarising or plotting contig information

    Depending on how the FastA file was created, different types of plots can be
    are available.  For instance, if the FastA was created with Canu,
    *nreads* and *covStat* information can be extracted. Therefore,
    plots such as :meth:`plot_scatter_contig_length_vs_nreads_cov`
    and :meth:`plot_contig_length_vs_nreads` can be used.

    """

    def __init__(self, filename, mode="canu"):

        """.. rubric:: **Constructor**

        :param filename: input FastA file
        :param canu: tool that created the output file.

        """
        super(Contigs, self).__init__(filename)
        self.mode = mode
        self._df = None

    def hist_plot_contig_length(self, bins=40, fontsize=16, lw=1):
        """Plot distribution of contig lengths

        :param bin: number of bins for the histogram
        :param fontsize: fontsize for xy-labels
        :param lw: width of bar contour edges
        :param ec: color of bar contours

        .. plot::

            from sequana import Contigs, sequana_data
            filename = sequana_data("test_contigs_spades.fasta")
            c = Contigs(filename)
            c.hist_plot_contig_length()

        """
        L = len(self.fasta.sequences)
        pylab.clf()
        pylab.hist(self.fasta.lengths, lw=lw, ec="k", bins=bins)
        pylab.grid()
        pylab.xlabel("Contig length", fontsize=fontsize)
        pylab.ylabel("#", fontsize=fontsize)
        pylab.title("Distribution {} contigs".format(L))

    def _get_df(self):
        if self._df is None:
            try:
                self._compute_spades_df()
            except ValueError:
                self._compute_df()
        return self._df

    df = property(_get_df)

    def _compute_spades_df(self):
        lengths = []
        names = []
        covs = []
        for name in self.fasta.names:
            _, ID, _, length, _, cov = name.split("_")
            lengths.append(length)
            names.append(ID)
            covs.append(cov)
        self._df = pd.DataFrame({"cov": covs, "length": lengths, "name": names})
        self._df = self._df.astype({"length": int, "cov": float})
        self._df = self._df[["name", "length", "cov"]]
        self._df["GC"] = self.get_gc()

    def _compute_df(self, window=100):
        data = tools._base_content(self.filename, window, "GC")
        names = self.fasta.names
        lengths = self.fasta.lengths
        GC = [np.nanmean(data[name]) for name in names]
        nreads = [0] * len(GC)
        covStats = [0] * len(GC)
        if self.mode == "canu":
            for i, comment in enumerate(self.fasta.comments):
                read = [x for x in comment.split() if x.startswith("reads")][0]
                covStat = [x for x in comment.split() if x.startswith("covStat")][0]
                read = read.split("=")[1]
                covStat = covStat.split("=")[1]
                nreads[i] = int(read)
                covStats[i] = float(covStat)
        df = pd.DataFrame({"GC": list(GC), "length": lengths, "name": names, "nread": nreads, "covStat": covStats})
        self._df = df.copy()
        return df

    def plot_contig_length_vs_nreads(self, fontsize=16, min_length=5000, min_nread=10, grid=True, logx=True, logy=True):
        """Plot contig length versus nreads

        In canu, contigs have the number of reads that support them.
        Here, we can see whether contigs have lots of reads supported them or not.

        .. note:: For Canu output only

        .. plot::

            from sequana import Contigs, sequana_data
            filename = sequana_data("test_contigs_ex1.fasta")
            c = Contigs(filename)
            c.plot_contig_length_vs_nreads(logx=False)

        """
        # same as plot_scatter_contig_length_nread_cov but no covStats information
        if not "nread" in self.df.columns:
            logger.warning("plot_scatter_contig_length_nread_cov required 'nread' column information (Canu output)")
            return
        pylab.clf()

        m1 = self.df.length.min()
        M1 = self.df.length.max()
        pylab.plot(self.df.length, self.df.nread, "o")
        pylab.xlabel("Contig length", fontsize=fontsize)
        pylab.ylabel("Contig N reads", fontsize=fontsize)
        if grid:
            pylab.grid()
        if logx:
            pylab.semilogx()
        if logy:
            pylab.semilogy

        query = "nread>@min_nread and length>@min_length"
        X = self.df.query(query)["length"]
        Y = self.df.query(query)["nread"]

        try:  # pragma: no cover
            A = np.vstack([X, np.ones(len(X))]).T
            m, c = np.linalg.lstsq(A, Y.as_matrix())[0]
            x = np.array([m1, M1])
            pylab.plot(x, m * x + c, "o-r")
        except AttributeError:
            pass

        pylab.xlabel("Contig length", fontsize=16)
        pylab.ylabel("nread support", fontsize=16)
        pylab.tight_layout()

    def plot_scatter_contig_length_vs_nreads_cov(
        self, fontsize=16, vmin=0, vmax=50, min_nreads=20, min_length=5000, grid=True, logx=True, logy=True
    ):
        """Scatter plot showing number of support reads and contig lengths

        .. note:: only for Canu output.

        .. plot::

            from sequana import Contigs, sequana_data
            filename = sequana_data("test_contigs_ex1.fasta")
            c = Contigs(filename)
            c.plot_scatter_contig_length_vs_nreads_cov()
        """
        if not "covStat" in self.df.columns:
            logger.warning(
                "plot_scatter_contig_length_nread_cov required 'covStat' coverage column information (Canu output). You may use plot_contig_length_vs_nreads method instead"
            )
            return

        if not "nread" in self.df.columns:  # pragma: no cover
            logger.warning("plot_scatter_contig_length_nread_cov required 'nread' column information (Canu output)")
            return

        m1 = self.df.length.min()
        M1 = self.df.length.max()

        # selection
        query = "nread>@min_nreads and length>@min_length"
        X = self.df.query(query)["length"]
        Y = self.df.query(query)["nread"]
        Z = self.df.query(query)["covStat"]

        if len(X) == 0:
            logger.warning("No contig after filtering. Set min_reads and min_length")
            return

        pylab.clf()
        pylab.scatter(X, Y, c=Z, vmin=vmin, vmax=vmax)
        pylab.colorbar()
        pylab.xlabel("Contig length", fontsize=fontsize)
        pylab.ylabel("Contig reads", fontsize=fontsize)

        try:  # pragma: no cover
            A = np.vstack([X, np.ones(len(X))]).T
            m, c = np.linalg.lstsq(A, Y.as_matrix())[0]
            x = np.array([m1, M1])
            pylab.plot(x, m * x + c, "o-r")
        except AttributeError:
            pass

        if grid:
            pylab.grid()
        if logx:
            pylab.semilogx()
        if logy:
            pylab.semilogy()

        pylab.tight_layout()

    """def get_contig_per_chromosome(self):
        if self.bam is None:
            print("no bam file found")
            return
        self.df = self.bam.get_df()
        df = self.df.query("flag in [0,16]")
        alldata = {}
        for chrom in sorted(df.rname.unique()):
            data = df.query("rname == @chrom").sort_values(by="rstart")[["qname", "qlen", "rstart", "rend"]]
            alldata[chrom] = data
        return alldata
    """

    def stats(self):
        """Return N50, L50 and total cumulated length"""
        from sequana.stats import N50, L50

        length = self.df["length"]
        return {"N50": N50(length), "total_length": sum(length), "L50": L50(length)}
