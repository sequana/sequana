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
"Pacbio amplicon related tools"
import glob

from sequana.freebayes_vcf_filter import VCF_freebayes, Variant
from sequana import BAM

import colorlog
logger = colorlog.getLogger(__name__)



import pandas as pd
import pylab



class Consensus():

    def __init__(self, bases, freebayes=None):
        self.filename_bases = bases
        self.min_depth = 10
        self.min_score = 1

        if freebayes is not None:
            v = VCF_freebayes(freebayes)
            self.variants = [Variant(x) for x in v if x]
        else:
            self.variants = []

    def identify_deletions(self):

        deletions = []
        for variant in self.variants:
            alt = variant.resume["alternative"]
            ref = variant.resume["reference"]
            pos = variant.resume['position']
            dp  = variant.resume['depth']
            score = variant.resume['freebayes_score']

            # filter out bad scores and bad depth
            if score < self.min_score:
                continue
            if dp < self.min_depth:
                continue

            # we need at least 2 bases for an insertion
            if len(alt) <1:
                continue

            # we need an alternate > than reference for an insertion
            if len(alt)<=len(ref):
                continue

            deletions.append(variant)
        return deletions

    def get_bases(self):

        # header of the consensus created by IGV, may have a warning in the
        # header, which should be ignored.
        toskip = 0
        with open(self.filename_bases) as fin:
            for line in fin.readline():
                if line.startswith(tuple(x for x in "1234567890")):
                    break
                else:
                    toskip += 1

        df = pd.read_csv(self.filename_bases, sep="\t", skiprows=toskip, header=None)
        df.columns = ["Pos", "A", "C", "G", "T", "N", "DEL", "INS"]
        df = df.set_index("Pos")

        # Low coverage values should be taken care of.
        # If below min_depth, we set everything to zero and consider the base as
        # an N. 
        indices = df[df.sum(axis=1)<self.min_depth].index
        df.loc[indices] = 0
        df.loc[indices, "N"] = 10000 # not important since we normalise later
        return df

    def run(self):

        # To normalise one need to ignore the insertions since there
        # are already included in the ACGT nucleotides
        cols = ["A", "C", "G", "T", "N", "DEL"] 
        df = self.get_bases()
        deletions = self.identify_deletions()

        # consensus without deletions
        dd = df.apply(lambda x: x.idxmax(), axis=1)

        # check that deletions are consistent with the data
        for d in deletions:
            pos = int(d.resume["position"])
            ref = d.resume["reference"]
            # compare the reference of the deletions with the consensus
            if "".join(dd.loc[pos:pos+len(ref)-1]) != ref:
                logger.warning("reference string {} not found in consensus at position {}".format(ref, pos))

        # Now, we insert the deletions removing the reference and then adding
        # the alternate
        # We aware that some deletions may overlap 
        for d in deletions:
            pos = int(d.resume["position"])
            ref = d.resume["reference"]
            alt = d.resume["alternative"]

            # the data up to the position of the reference/alternate SNP
            # indices may not start at zero so we use loc instead of iloc
            dfA = df.loc[0:pos-1]

            # the alternate data needs a dummy dataframe. The indices are 
            # e.g. 0,1,2,3,4,5 and 
            # We reset the indices to start at the dfA last position and
            # to be constant. For instance a dataframe dfB of 3 rows to be
            # appended aftre position 1500 will have the indices 1500,1500,1500
            # This garantee that the following dataframe dfC has indices > to
            # those of dfB while allowing the next iteration to use the same
            # consistance indices when searching for the next deletions
            dfB = df.iloc[0:len(alt)].copy()
            dfB.index = [pos] *  len(dfB)
            dfB *= 0
            for i, nucleotide in enumerate(alt):
                dfB.iloc[i][nucleotide] = 10000 

            # the rest of the data
            dfC = df.loc[pos+len(ref):]

            # !! do no reset indices !!! so that inserted dfB is still sorted
            # and next accession with iloc/loc are still correctin the next
            # iteration
            df = dfA.append(dfB).append(dfC)#.reset_index(drop = True)

        # now we can reset the indices
        df.reset_index(drop=True, inplace=True)

        dd = df.apply(lambda x: x.idxmax(), axis=1)
        return dd

    def save_consensus(self, output, identifier):

        dd = self.run()
        with open(output, "w") as fout:
             data = "".join(dd).replace("DEL", "")
             sample = identifier
             fout.write(">{}\n{}\n".format(sample,data))

    def get_population(self, threshold=0.1, Npop=2):
        df = self.get_bases()
        cols = ["A", "C", "G", "T", "N", "DEL"]
        df = df.divide(df[cols].sum(axis=1), axis=0)
        selection = df[(df[cols]>threshold).sum(axis=1)>=2]
        return selection




class LAA():
    """Reads a set of file called amplicon_summary.csv generated by the LAA
    pipeline from smrtlink pacbio tool.

    The file have this format::

        BarcodeName,FastaName,CoarseCluster,Phase,TotalCoverage,SequenceLength,PredictedAccuracy,ConsensusConverged,NoiseSequence,IsDuplicate,DuplicateOf,IsChimera,ChimeraScore,ParentSequenceA,ParentSequenceB,CrossoverPosition
        24--24,Barcode24--24_Cluster0_Phase0_NumReads497,0,0,497,2958,1,1,0,0,N/A,0,-1,N/A,N/A,-1
        24--24,Barcode24--24_Cluster1_Phase0_NumReads496,1,0,496,3137,1,1,0,0,N/A,0,-1,N/A,N/A,-1

    See :meth:`plot_max_length_amplicon_per_barcode` for some details about
    sample names.

    """
    def __init__(self, where="bc*"):
        self.filenames = glob.glob(where + "/" + "amplicon_*summary.csv")
        self.data = [pd.read_csv(this) for this in self.filenames]
        self.numbers = [len(x) for x in self.data]

    def hist_amplicon(self, fontsize=12):
        pylab.hist(self.numbers, bins=max(self.numbers), ec="k", align="left")
        pylab.ylabel("#", fontsize=fontsize)
        pylab.ylabel("Number of amplicons per barcode", fontsize=fontsize)

    def plot_max_length_amplicon_per_barcode(self, sample_names=None):
        """Plot max length of the amplicons per barcode

        :param sample_names: names of the barcode. If not provided use number
            from 1 to N. See below.

        One difficulty here is to associate a file with a barcode name. The files
        created by LAA have the same basename and the content of the file cannot be
        parsed to obtain the barcode (indeed some may be empty). So, we need to
        provide a list of sample names in :meth:`plot_max_length_amplicon_per_barcode` 
        associated to the files. Otherwise a simple range of names from 1 to N is
        used.
        """
        if sample_names is None:
            sample_names = range(1, len(self.data)+1)

        data = [max(x.SequenceLength) if len(x.SequenceLength) else 0  for x in self.data]
        pylab.plot(sample_names, data, "o")
        pylab.xlabel("barcode name (filename)")
        pylab.ylabel("Max length amongst amplicons (0 if no amplicon)")
        #savefig("max_amplicon_length_per_barcode.png", dpi=200)


class LAA_Assembly():
    """

    Input is a SAM/BAM from the mapping of amplicon onto a known reference.
    Based on the position, we can construct the new reference.

    """
    def __init__(self, filename):
        self.bam = BAM(filename)


    def build_reference(self):
        self.bam.reset()
        # scan BAM file assuming it is small
        aa = [a for a in self.bam]

        # retrieve data of interest
        data = [(a.pos, {
                    "name":a.query_name,
                    "sequence": a.query_sequence,
                    "cigar": a.cigarstring,
                    "position": a.pos,
                    "qstart": a.qstart,
                    "qend": a.qend}) for a in aa]

        # sort by starting position
        data.sort(key=lambda x: x[0])

        for i, read in enumerate(data):
            read = read[1]
            if i == 0:
                sequence = read["sequence"]     # 2 is query_sequence
            else:
                pr = data[i-1][1]   # previous read
                L = len(pr["sequence"])
                end_position_pr = pr['position'] - pr['qstart'] + L 

                # overlap between previous read and this one
                overlap = end_position_pr - (read['position'] - read['qstart']) +0
                print(overlap)
                print(pr['position'], pr['qstart'], L, end_position_pr)
                print(read['position'], read['qstart'])
                sequence = sequence + read["sequence"][overlap+1:]

        # argmax([sum(a==b for a,b in zip(X[-i:] , Y[:i]))/float(i+1) for i in range(1000)])
        return sequence

    def save_fasta(self, filename, sequence=None):
        if sequence is None:
            sequence = self.build_reference()


        with open(filename, "w") as fout:
            fout.write(">test\n{}".format(sequence))


