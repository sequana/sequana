from sequana.lazy import pylab
from sequana import FastA, BAM
from sequana.lazy import pandas as pd
from sequana.lazy import numpy as np
from sequana import tools


class Contigs(object):

    def __init__(self, filename, reference=None, bamfile=None, mode="canu"):

        """


            minimap2 -x map-pb reference filename -a > temp.sam
            bioconvert sam2bam temp.sam temp.bam

        """
        self.filename = filename
        self.fasta = FastA(filename)
        self.mode = mode
        self._df = None
        if bamfile:
            self.bam = BAM(bamfile)
        else:
            self.bam = None
        self.reference = reference

    def bar_plot_contigs_length(self):
        # show length of N contigs as compare to length of the reference
        fref = FastA(self.reference)
        Nref = len(fref.sequences)
        N = len(self.fasta)
        pylab.clf()
        pylab.bar(range(0, N, int(pylab.ceil(N/Nref))), sorted(fref.lengths), width=Nref/1.1,
            label="Plasmodium chromosomes")
        pylab.bar(range(0, N), sorted(self.fasta.lengths), width=1,
            label="canu {} contigs".format(N))
        pylab.legend()
        #pylab.savefig("1179_195_contigs.png", dpi=200)

    def hist_plot_contig_length(self, bins=40, fontsize=16):
        """Plot distribution of contig lengths"""
        L = len(self.fasta.sequences)
        pylab.hist(self.fasta.lengths, lw=1, ec="k",bins=bins) 
        pylab.grid()
        pylab.xlabel("Contig length", fontsize=fontsize)
        pylab.ylabel("#", fontsize=fontsize)
        pylab.title("Distribution {} contigs".format(L))

    def get_df(self, window=100):
        print("building GC content")
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
        #if self.bamfile
        df = pd.DataFrame({"GC":list(GC), "length":lengths, "name": names,
                           "nread": nreads, "covStat": covStats })

        # deal with the bamfile
        if self.bam:
            bam_df = self.bam.get_df()
            bam_df = bam_df.query("flag in [0,16]")
            bam_df.set_index("qname", inplace=True)
            chrom_name = bam_df.loc[self.fasta.names]["rname"]
            df["chromosome"] = list(chrom_name)

        self._df = df.copy()
        return df

    def plot_contig_length_vs_nreads(self, fontsize=16):
        # same as plot_scatter_contig_length_nread_cov
        if self._df is None:
            _ = self.get_df()
        pylab.clf()
        df = self._df

        m1 = df.length.min()
        M1 = df.length.max()
        pylab.loglog(df.length, df.nread, "o")
        pylab.xlabel("Contig length", fontsize=fontsize)
        pylab.ylabel("Contig N reads", fontsize=fontsize)
        pylab.grid()

        X = df.query("nread>10 and length>100000")['length']
        Y = df.query("nread>10 and length>100000")['nread']
        A = np.vstack([X, np.ones(len(X))]).T
        m, c = np.linalg.lstsq(A, Y.as_matrix())[0]
        x = np.array([m1, M1])
        pylab.plot(x, m*x+c, "o-r")
        pylab.tight_layout()

    def plot_scatter_contig_length_nread_cov(self, fontsize=16, vmin=0, vmax=50,
            min_nreads=20, min_length=50000):

        if self._df is None:
            _ = self.get_df()
        pylab.clf()
        df = self._df

        m1 = df.length.min()
        M1 = df.length.max()

        # least square
        X = df.query("nread>@min_nreads and length>@min_length")['length']
        Y = df.query("nread>@min_nreads and length>@min_length")['nread']
        Z = df.query("nread>@min_nreads and length>@min_length")['covStat']
        print(X)
        print(Y)
        print(Z)

        A = np.vstack([X, np.ones(len(X))]).T
        m, c = np.linalg.lstsq(A, Y.as_matrix())[0]
        x = np.array([m1, M1])

        X = df['length']
        Y = df['nread']
        Z = df['covStat']
        pylab.scatter(X, Y, c=Z, vmin=vmin, vmax=vmax)
        pylab.colorbar()
        pylab.xlabel("Contig length", fontsize=fontsize)
        pylab.ylabel("Contig reads", fontsize=fontsize)
        pylab.title("coverage function of contig length and reads used")
        pylab.grid()
        pylab.plot(x, m*x+c, "o-r")
        pylab.loglog()
        pylab.tight_layout()

    def get_contig_per_chromosome(self):
        if self.bam is None:
            print("no bam file found")
            return
        df = self.bam.get_df()
        df = df.query("flag in [0,16]")
        alldata = {}
        for chrom in sorted(df.rname.unique()):
            data = df.query("rname == @chrom").sort_values(by='rstart')[["qname", "qlen", "rstart","rend"]]
            alldata[chrom] = data
        return alldata

    def stats(self):
        from sequana.stats import N50
        length = self.get_df()['length']
        return {
            'N50': N50(length),
            'total_length': sum(length),
            'L50': L50(length)
            }

    def plot_contig_length_vs_GC(self):
        pylab.plot(self.get_df()["length"], self.get_df()['GC'], "o")

