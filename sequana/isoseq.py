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
import glob
import os
from collections import Counter


from sequana.lazy import  pandas as pd
from sequana.lazy import  numpy as np
from sequana import FastQ, FastA
from sequana.pacbio import PacbioSubreads
from sequana.lazy import pylab
from sequana import phred


import colorlog
logger = colorlog.getLogger(__name__)



__all__ = ["IsoSeqQC", "IsoSeqBAM"]


class IsoSeqQC(object):
    """


    Use get_isoseq_files on smrtlink to get the proper files

    iso = IsoSeqQC()
    iso.hist_read_length_consensus_isoform() # histo CCS
    iso.stats() # "CCS" key is equivalent to summary metrics in CCS report


    todo: get CCS passes histogram . Where to get the info of passes ? 

    """
    def __init__(self, directory=".", prefix=""):
        self.prefix = prefix
        self.directory = directory
        self.sample_name = "undefined"

        # low quality isoforms
        filename = "all.polished_lq.fastq"
        self.lq_isoforms = self.get_file(filename)
        if self.lq_isoforms:
            logger.info("Reading {}".format(filename))
            self.lq_sequence = FastQ(self.lq_isoforms)

        # high quality isoforms
        filename = "all.polished_hq.fastq"
        self.hq_isoforms = self.get_file(filename)
        if self.hq_isoforms:
            logger.info("Reading {}".format(filename))
            self.hq_sequence = FastQ(self.hq_isoforms)

        # General info
        filename = "file.csv"
        self.csv = self.get_file(filename)
        if self.csv:
            logger.info("Reading {}".format(filename))
            self.data = pd.read_csv(self.csv)

        # CCS fasta sequence
        #self.ccs = self.get_file("-ccs.tar.gz")
        filename = "ccs.fasta"
        self.ccs = self.get_file(filename, noprefix=True)
        if self.ccs:
            logger.info("Reading {}".format(filename))
            self.ccs = FastA(self.ccs)

    def get_file(self, tag, noprefix=False):
        if noprefix:
            filenames = glob.glob(self.directory + os.sep + tag)
        else:
            filenames = glob.glob(self.directory + os.sep + self.prefix + tag)
        if len(filenames) == 1:
            return filenames[0]
        elif len(filenames) > 1:
            print("Found several files ending in %s" % tag)
        else:
            print("No files matching %s" % tag)
        return None

    def stats(self):
        results = {}
        if self.data is not None:
            logger.info("Reading strand")
            results['strand'] = {
                "+": sum(self.data.strand == "+"),
                "-": sum(self.data.strand == "-"),
                "?": sum(self.data.strand.isnull())
            }

            results['classification'] = {
                "total_ccs_reads" : len(self.data),
                "five_prime_reads" : int(self.data.fiveseen.sum()),
                "three_prime_reads" : int(self.data.threeseen.sum()),
                "chimera" : int(self.data.chimera.sum()),
                "polyA_reads" : int(self.data.polyAseen.sum()),
            }

        if self.lq_isoforms:
            logger.info("Reading LQ isoforms")
            results['lq_isoform'] = self.lq_sequence.stats() # number of 

        if self.hq_isoforms:
            logger.info("Reading HQ isoforms")
            results['hq_isoform'] = self.hq_sequence.stats() # number of polished HQ isoform

        if self.ccs:
            seq = [ len(read.sequence) for read in self.ccs]
            results["CCS"] = {
                "mean_length" : pylab.mean(seq),
                "number_ccs_bases" : sum(seq),
                "number_ccs_reads" : len(seq)
            }

        self.idents_v = []
        self.full_v = []
        self.non_full_v = []
        self.isoform_lengths = []
        for read in self.lq_sequence:
            ident, full, non_full, length = read['identifier'].decode().split(";")
            self.idents_v.append(ident)
            self.full_v.append(int(full.split("=")[1]))
            self.non_full_v.append(int(non_full.split("=")[1]))
            self.isoform_lengths.append(int(length.split("=")[1]))

        return results

    def to_summary(self, filename="sequana_summary_isoseq.json", data=None):
        """Save statistics into a JSON file

        :param filename:
        :param data: dictionary to save. If not provided, use :meth:`stats`

        """
        from sequana.summary import Summary
        if data is None:
            data = self.stats()
        Summary("isoseq",self.sample_name, data=data).to_json(filename)

    def hist_read_length_consensus_isoform(self, mode="all", bins=80, rwidth=0.8,
        align="left", fontsize=16, edgecolor="k", **kwargs):
        """

        mode can be all, lq, hq
        """
        pylab.clf()

        L1 = [len(read['sequence']) for read in self.lq_sequence]
        L2 = [len(read['sequence']) for read in self.hq_sequence]
        if mode == "all":
            L = L1 + L2
        elif mode == "lq":
            L = L1
        else:
            L = L2
 
        Y, X, _ = pylab.hist(L, bins=bins, rwidth=rwidth, align=align,
                    ec=edgecolor, **kwargs)
        pylab.gca().set_ylim(bottom=0)
        pylab.gca().set_xlim(left=0)
        pylab.xlabel("Read length", fontsize=fontsize)
        pylab.ylabel("Number of reads", fontsize=fontsize)
        pylab.grid()

        ax_twin = pylab.gca().twinx()

        shift = (X[1] - X[0]) / 2

        ax_twin.plot(X[0:-1]- shift, len(L) - pylab.cumsum(Y), "k")
        ax_twin.set_ylim(bottom=0)
        pylab.ylabel("CDF", fontsize=fontsize)

        pylab.title("Read length of Consensus isoforms reads")

    def hist_average_quality(self, fontsize=16, bins=None):
        """

        bins is from 0 to 94 
        """

        hq_qv = [pylab.mean([ord(X)-33 for X in read['quality'].decode()]) 
                for read in self.hq_sequence]
        lq_qv = [pylab.mean([ord(X) -33 for X in read['quality'].decode()]) 
            for read in self.lq_sequence]

        if bins is None:
            bins = range(0,94)
        Y1, X = np.histogram(hq_qv, bins=bins)
        Y2, X = np.histogram(lq_qv, bins=bins)
        pylab.bar(X[1:], Y1, width=1, label="HQ")
        pylab.bar(X[1:], Y2, bottom=Y1, width=1, label="LQ")
        pylab.xlim([0.5, 93.5])

        pylab.xlabel("Isoform average QV")
        pylab.ylabel("# Isoform")
        pylab.legend(fontsize=fontsize)

        ax = pylab.twinx()
        N = np.sum(Y1+Y2)
        ax.plot(X, [N] + list(N-np.cumsum(Y1+Y2)), "k")


class IsoSeqBAM(object):
    """Reads raw IsoSeq BAM file (subreads)

    Creates some plots and stats

    """
    def __init__(self, filename):
        self.filename = filename
        self._passes = None
        self._lengths = None

    @property
    def lengths(self):
        if self._lengths is None:
            _ = self.passes # extract passes and lengths
        return self._lengths

    @property
    def passes(self):
        from collections import Counter
        if self._passes is None:
            # for BAM of 1M reads, takes 10-15 seconds
            from pysam import AlignmentFile
            ccs = AlignmentFile(self.filename, check_sq=False)
            qnames = [a.qname for a in ccs]
            names = [int(qname.split("/")[1]) for qname in qnames]
            self.qnames = qnames
            lengths = []
            for qname in qnames:
                a,b = qname.split("/")[2].split("_")
                lengths.append(int(b)-int(a))
            self._lengths = lengths
            self._passes = list(Counter(names).values())
        return self._passes

    def hist_passes(self, bins=100):
        pylab.hist(self.passes, bins=bins)

    def hist_read_length(self, bins=100):
        pylab.hist(self.lengths, bins=bins)

    def stats(self):
        return {"mean_read_length": pylab.mean(self.lengths), 
                "ZMW": len(self.passes),
                "N": len(self.lengths),
                "mean_ZMW_passes": pylab.mean(self.passes)}


class PacbioIsoSeqMappedIsoforms(object):
    """Here, we load a SAM/BAM file generated with minimap using
    as input the BAM file  created with the mapping og HQ isoforms
    on a reference.

    df contains a dataframe for each read found in the SAM (and hq_isoform)
    we populate the GC content, the mapping flag, the reference name (-1 means
    no mapping i.e flag ==4). flag of 4 means unmapped and there is no
    ambiguity about it.

    In the data file example, other flags are 0, 16 (SEQ being reverse
    complement<F12>) , 2048 (supplementary     segment).

    Example of minimap2 command::

        minimap2 -t 4  -ax splice -uf --secondary=no  SIRV-E0.fa
            hq_isoforms.fasta 1> hq_isoforms.fasta.sam 2> hq_isoforms.fasta.sam.log

    Reads a SAM file for now. BAM should work as well

    """
    def __init__(self, filename):
        self.filename = filename
        self.bam = PacbioSubreads(self.filename)
        self._df = None

    @property
    def df(self):
        if self._df is not None:
            return self._df

        # !! for isoseq, we should be able to load everything into memory
        self.bam.reset()
        data = [a for a in self.bam.data]

        df = self.bam.df.copy()

        rnames = [self.bam.data.get_reference_name(a.rname) if a.rname!=-1
                  else -1 for a in data]

        df['reference_name'] = rnames
        df['flags'] = [a.flag for a in data]
        df['mapq'] = [a.mapq for a in data]
        df['cigar'] = [a.cigarstring for a in data]
        df['qname'] = [a.qname for a in data]

        # Drop SNR that are not populated in the mapped BAM file.
        df.drop(['snr_A', 'snr_C', 'snr_G', 'snr_T'], axis=1, inplace=True)


        # TODO. input could be basde on mapping of CCS in which case, the ZMW is
        # stored and the following does not work. could check whether the
        # pattern is pXXfXX
        try:
            df["full_length"] = df["qname"].apply(lambda x: int(x.split('/')[1].split("p")[0].strip("f")))
            df["non_full_length"] = df["qname"].apply(lambda x: int(x.split("/")[1].split("p")[1].strip("f")))
        except:
            pass

        self._df = df

        return self._df

    def hist_isoform_length_mapped_vs_unmapped(self, bins=None):
        df = self.df
        if bins is None:
            bins = range(0, len(df.reference_length.max()), 100)
        mapped = df[df.reference_name != -1]
        unmapped = df[df.reference_name == -1]
        pylab.hist(mapped.reference_length, bins=bins, alpha=0.5,
            label="mapped {}".format(len(mapped)), density=False)
        pylab.hist(unmapped.reference, bins=bins, alpha=0.5,
            label="unmapped {}".format(len(unmapped)), density=False)
        pylab.xlabel("Isoform length")
        pylab.legend()

    def hist_transcript(self, hide_unmapped=True):
        pylab.clf()

        if hide_unmapped is True:
            query = "reference_length>0 and reference_name!=-1"
        else:
            query = "reference_length>0"

        print(query)
        ts = self.df.query(query).groupby("reference_name").count().reference_length
        if len(ts) == 0:
            print("nothing to plot")
            return ts

        ts.plot(kind="bar" ,color="r")
        try: pylab.tight_layout()
        except: pass
        return ts

    def bar_mapq(self, logy=True, xmin=0, xmax=60, fontsize=12):
        self.df.mapq.hist()
        if logy:
            pylab.semilogy()
        pylab.xlim([xmin, xmax])
        pylab.xlabel("Mapping quality", fontsize=fontsize)

    def plot_sirv_by_group(self, title, shift=5, plot=False, mapq_min=-1):
        aa = self.df.query("reference_name not in [-1, '-1']").copy()
        if len(aa) == 0:
            return pd.Series(), self.df

        aa['group'] = aa.reference_name.apply(lambda x: x[0:shift])
        mapped = aa.query("mapq>@mapq_min").groupby("group").count()["mapq"]
        mapped.name = None

        if plot:
            mapped.plot(kind="bar")
            pylab.title(title)
            pylab.tight_layout()
        #data.to_csv(path + "_hq_sirv_grouped.csv")
        return mapped, self.df


class SIRVReference():
    """

    See https://www.lexogen.com/sirvs/downloads/ to download one of the XLSX
    file

    ::

        ss = isoseq.SirvReference()
        ss.from_excel("SIRV_Set2.xls")
        ss.to_fasta("temp.fasta")


    """
    def __init__(self):
        pass

    def from_excel(self, filename, nrows=100):
        """
        skip 4 empty rows. next two are a header split on two lines

        """
        # header is split on two rows, so we need to skip it 
        df = pd.read_excel(filename,  skiprows=4, header=None)
        assert df.iloc[1, 4] == "c" # this column will be used
        assert df.iloc[0, 2] == "SIRV ID"
        assert df.iloc[0, 1] == "SIRV gene"
        assert df.iloc[1, 23] == "length (nt)"
        assert df.iloc[1, 25] == "pA length (nt)"
        assert df.iloc[1, 26] == "Sequence" # with polyA
        assert df.iloc[1, 29] == "Sequence" # without polyA !! 
        # in excel this last column is actually an equation based on the
        # sequence with polyA and may be empty when using from_excel method.

        # filter out useless columns and first two lines
        df = df.iloc[2:, [2,4,23,25,26, 29]]
        # drop useless rows !! keep this statement after the first filtering
        df = df[df.iloc[:,0].isnull() == False]

        # drop header
        df.reset_index(inplace=True, drop=True)


        df.columns = ["sirv_id", "annotation", "length",
                      "length_polyA", "sequence_with_polya", "sequence"]

        # remove and replace sequence column using sequence_with_poly and
        # removing the polyA
        #df.drop("sequence", inplace=True)
        seq = [a[0:-b] for a,b in zip(df.sequence_with_polya, df.length_polyA)]
        df['sequence'] = seq


        self.df = df

    def to_fasta(self, filename):
        data = self.df.query("annotation == 1").copy()
        with open(filename, "w") as fout:
            for _, this in data.iterrows():
                fout.write(">{}\n{}\n".format(this.sirv_id, this.sequence))


class SIRV(object):
    def __init__(self, filename, shift=4):
        from sequana import FastA
        self.SIRV = FastA(filename)
        self.shift = 5

    def __len__(self):
        return len(self.SIRV.names)

    @property
    def group_names(self):
        data = list(set([x[0:self.shift] for x in self.SIRV.names]))
        return sorted(data)

    @property
    def group_lengths(self):
        data = dict(Counter( list([x[0:self.shift] for x in self.SIRV.names])))
        return data


class PacbioIsoSeqMultipleIsoforms(object):
    """

        mul.read("./180223_113019_lexo1/data/ccs_sirv.sam", "lexo1 ccs")
        mul.read("./180223_221457_lexo2/data/ccs_sirv.sam", "lexo2 ccs")
        mul.read("./180206_193114_lexo3/data/ccs_sirv.sam", "lexo3 ccs")
        mul.read("./180224_083446_lexo3/data/ccs_sirv.sam", "lexo3 ccs bis")

    """
    SIRV_names = ["SIRV1", "SIRV2", "SIRV3", "SIRV4", "SIRV5", "SIRV6", "SIRV7"]
    SIRV_total_lengths = [7075, 6395, 12988, 5592, 16206, 15278, 13375]
    def __init__(self, sirv=None, sirv_shift=5):
        self.labels = []
        self.rawdata = []
        if sirv:
            try:
                self.SIRV_data = SIRV(sirv, shift=sirv_shift)
                self.SIRV_names = self.SIRV_data.group_names
                self.SIRV_lengths = self.SIRV_data.SIRV.get_lengths_as_dict()
            except:
                print("Could not read {}".format(sirv))

        self.mapq_min = 0

    @property
    def N(self):
        return [len(x) for x in self.rawdata]

    @property
    def sirv(self):
        sirv = []
        shift=5
        for df in self.rawdata:
            aa = df.query("reference_name not in [-1, '-1']").copy()
            if len(aa) == 0:
                sirv.append(pd.Series(index=self.SIRV_names))
                continue

            aa['group'] = aa.reference_name.apply(lambda x: x[0:shift])
            # filter quality and flags
            #mask = np.logical_or(df.flags & 256, df.flags & 2048)
            #aa = aa.query("mapq>=@self.mapq_min and @mask")

            data = aa.groupby("group").count()["mapq"]
            data.name = None
            sirv.append(data)
        return sirv

    def filter_raw_data(self, quality=0, flags=[256,2048]):
        for i, df in enumerate(self.rawdata):
            mask = np.logical_or(df['flags'] & 256, df['flags'] & 2048)
            ones = [-1, '-1']
            df = df.query("mapq>=@quality and not @mask and reference_name not in @ones")
            self.rawdata[i] = df

    def read(self, filename, tag):
        if filename.endswith(".csv"):
            data = pd.read_csv(filename)
            data = data.query("reference_name not in [-1, '-1']")
        else:
            sam = PacbioIsoSeqMappedIsoforms(filename)
            mapped, data = sam.plot_sirv_by_group(tag)

        self.rawdata.append(data)
        self.labels.append(tag)

    def plot_corr(self):
        lengths = self.SIRV_data.SIRV.get_lengths_as_dict()
        spikes = self.spikes_found()
        spikes["lengths"] = [lengths[k] for k in spikes.index]
        corr = spikes.corr()
        pylab.imshow(corr)
        N = len(spikes.columns)
        pylab.xticks(range(N), spikes.columns, rotation=90)
        pylab.yticks(range(N), spikes.columns)
        pylab.clim(0,1)
        pylab.colorbar()

    def plot_bar_grouped(self, normalise=False, ncol=2, N=None):
        """

        :param normalise:
        :param ncol: columns in the legend

        """
        if N is not None:
            N = np.array(N)
        else:
            N = np.array([len(x) for x in self.rawdata])

        dd = pd.DataFrame(self.sirv).T
        if normalise:
            dd = dd/ (N/max(N))
        dd.columns = self.labels

        dd.plot(kind="bar")
        pylab.xlabel("")
        pylab.legend(self.labels, ncol=ncol)
        pylab.tight_layout()
        return dd

    def spikes_found(self, spikes_filename=None):
        if spikes_filename:
            sirv_names = SIRV(spikes_filename).SIRV.names
        else:
            sirv_names = []
            for data in self.rawdata:
                sirv_names.extend(data.reference_name.unique())
            sirv_names = set(sirv_names)
            try:sirv_names.remove(-1)
            except:pass
            sirv_names = sorted(list(sirv_names))

        df = pd.DataFrame(index=sirv_names, columns=self.labels)
        for i, data in enumerate(self.rawdata):
            aa = data.query("reference_name not in [-1, '-1']").copy()
            sirv_detected = aa.groupby("reference_name").count()['mapq']
            df[self.labels[i]] = sirv_detected
        return df

    def plot_bar(self, spikes_filename=None, ratio=100):
        data = self.spikes_found(spikes_filename)
        lengths = [self.SIRV_lengths[x] for x in data.index]
        data.plot(kind="bar")
        pylab.plot(np.array(lengths)/ratio)
        pylab.tight_layout()
        return data




def get_stats(pattern, mode):
    data = [[],[]]
    for filename in glob.glob(pattern + "/data/{}_sirv.sam".format(mode)):
        b = isoseq.PacbioIsoSeqMappedIsoforms(filename)
        N = Summary(filename.split("/")[0]+"/sequana_summary_isoseq.json").data["hq_isoform"]['N']
        data[0].append(N)
        data[1].append(len(b.df.query("reference_name not in [-1, '-1']")))
        print("scanned {}".format(filename))
    return data




class GeneCoverage():
    """

    - column 12 should be the percentage of gene coverage
    
    """
    def __init__(self, filename):
        self.filename = filename
        self.df = pd.read_csv(filename, header=None, sep="\t")
        self.coverage_column = 12

    def __str__(self):
        icol = self.coverage_column
        L = float(len(self.df))
        S0 = sum(self.df[icol]==0)
        S2 = sum(self.df[icol]==1)
        S1 = L - S0 - S2
        S90 = sum(self.df[icol]>0.9)
        S50 = sum(self.df[icol]>0.5)
        S99 = sum(self.df[icol]>0.99)
        txt = "Number of genes: {}\n".format(len(self.df))
        txt += "Fully detected genes: {} ({:.2f}%)\n".format(S2, S2/L*100)
        txt += "Partially detected genes: {} ({:.2f}%)\n".format(S1, S1/L*100)
        txt += "Detected genes (99% covered): {} ({:.2f}%)\n".format(S99, S99/L*100)
        txt += "Detected genes (90% covered): {} ({:.2f}%)\n".format(S90, S90/L*100)
        txt += "Detected genes (50% covered): {} ({:.2f}%)\n".format(S50, S50/L*100)
        txt += "Undetected genes: {} ({:.2f}%)\n".format(S0, S0/L*100)
        return txt

    def plot(self, X=[0,0.1,0.2,0.3,.4,.5,.6,.7,.8,.9,.95,.99,.999,1],
            fontsize=16, label=None):
        """plot percentage of genes covered (y axis) as a function of percentage
        of genes covered at least by X percent (x-axis). 

        """
        icol = self.coverage_column
        N = float(len(self.df))
        X = np.array(X) 
        Y = np.array([sum(self.df[icol]>x)/N*100 for x in X])
        if label is None:
            pylab.plot(X*100, Y, "o-")
        else:
            pylab.plot(X*100, Y, "o-", label=label)
        pylab.xlabel("Gene coverage (%)", fontsize=fontsize)
        pylab.ylabel("Percentage of genes covered", fontsize=fontsize)
        for this in [25,50,75]:
            pylab.axhline(this, color="r", alpha=0.5, ls="--")
            pylab.axvline(this, color="r", alpha=0.5, ls="--")

    def get_percentage_genes_covered_at_this_fraction(self, this):
        assert this<=1 and this>=0
        icol = self.coverage_column
        X = pylab.linspace(0, 1, 101)
        N = float(len(self.df))
        Y = np.array([sum(self.df[icol]>x)/N*100 for x in X])
        return np.interp(this, X, Y)








