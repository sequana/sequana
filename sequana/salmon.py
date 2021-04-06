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

import pandas as pd

from sequana.gff3 import GFF3

import colorlog
logger = colorlog.getLogger(__name__)


class Salmon():
    """Factory to read counts from salmon and create feature counts usable for deseq2"""

    def __init__(self, filename, gff_input, attribute="transcript_id"):
        self.filename = filename
        df = pd.read_csv(filename, sep='\t')

        self.df = df


        logger.info("Initialisation the GFF")
        # handling the GFF input file
        # The input gff may already be an instance off GFF
        if isinstance(gff_input, GFF3):
            self.gff = gff_input
        else:
            self.gff = GFF3(gff_input)
        # create the dataframe once for all
        self.gff.df

        # Of course, there is one gene by transcript, several transcript per gene
        # trs2genes is a dict with one transcript as key and one gene values
        # genes2trs is a dict with one gene sas key and a list of transcripts

    def get_feature_counts(self, feature=None, attribute=None):

        if self.df.Name.iloc[0].startswith("transcript:"):
            logger.info("salmon results start with transcript: tag. Eukaryotes mode")
            logger.info("Identifying the mapping transcript vs genes")
            self.trs2genes, self.genes2trs = self.gff.transcript_to_gene_mapping(attribute="transcript_id")
            results = self.get_feature_counts_eukaryotes(feature, attribute)
        else:
            logger.info("salmon results not starting with transcript. Prokaryotes mode")
            results = self.get_feature_counts_prokaryotes(feature, attribute)
        return results

    def get_feature_counts_eukaryotes(self, feature=None, attribute=None):

        if feature is None:
            feature = "gene"

        if attribute is None:
            attribute = "ID"

        # just to not loose the original
        df = self.df.copy()

        # Name contains the salmon entries read from gffread that uses
        # transcript_id. From this transcript id, we get the gene (parent)
        df['Gene'] = [self.trs2genes[x] for x in self.df.Name]

        #groups = df.groupby('Gene').groups
        counts_on_genes = df.groupby('Gene').NumReads.sum()


        ff = self.filename.split("/")[-1]
        results = f"\nGeneid\tChr\tStart\tEnd\tStrand\tLength\t{ff}"

        # mouse 25814 gene (feature)
        #       53715 gene_id (attribute)
        #      135181 transcript_id (attribute)
        #      133618 transcript_id from salmon
        #      135181 entries in transcript fasta (gffread)

        # gffread extact transcript_id from the gff if present
        # otherwise, extract geneID or gene_id
        logger.info("Recreating the feature counts")

        genes = {}

        dd = self.gff.df.query("ID in @counts_on_genes.index")
        dd = dd.set_index("ID")
        dd = dd.loc[counts_on_genes.index]
        self.dd = dd

        types = dd['type'].values
        starts = dd['start'].values
        stops = dd['stop'].values
        strands = dd['strand'].values
        seqids = dd['seqid'].values

        from easydev import Progress
        pb = Progress(len(counts_on_genes))

        S = 0

        logger.info("Grouping")
        TPMgroup = df.groupby('Gene').apply(lambda group: group['TPM'].sum())
        efflength_null = df.groupby('Gene').apply(lambda group: group['EffectiveLength'].mean())

        groups = df.groupby('Gene')
        for i, name in enumerate(counts_on_genes.index):
            # Since we use ID, there should be only one hit. we select the first
            # one to convert to a Series

            tpm_sum = TPMgroup.loc[name]
            if tpm_sum == 0:
                length = efflength_null.loc[name]
            else:
                abundances = groups.get_group(name).TPM
                efflength = groups.get_group(name).EffectiveLength
                length = sum([x*y for x,y in zip(abundances, efflength)]) / abundances.sum()
                S += abundances.sum()

            # FIXME we keep only types 'gene' to agree with output of
            # start/bowtie when working on the gene feature. What would happen
            # to compare salmon wit other type of features ? 
            if types[i] == "gene":
                start = starts[i]
                stop = stops[i]
                seqid = seqids[i]
                strand = strands[i]
                NumReads = counts_on_genes.loc[name]
                length = length
                name = name.replace("gene:", "")
                results += f"\n{name}\t{seqid}\t{start}\t{stop}\t{strand}\t{length}\t{NumReads}"
            else:
                pass
            pb.animate(i)
        return results
        """

In [179]: genes2trs['gene:ENSMUSG00000000028']
Out[179]: 
['transcript:ENSMUST00000000028',
 'transcript:ENSMUST00000096990',
 'transcript:ENSMUST00000115585']

length = 52660
sum(counts) = 24116676   OKAY if we include all features
sum(length) = 66698347   difference 66699167 of 820.1839 ???
sum(abundance) = 1e6

>head(counts)
ene:ENSMUSG00000000001 3266.000
gene:ENSMUSG00000000003    0.000
gene:ENSMUSG00000000028   68.000
gene:ENSMUSG00000000031 1113.000
gene:ENSMUSG00000000037  113.999
gene:ENSMUSG00000000049    1.000
> head(txi$length)   --> effective length average
                            [,1]
gene:ENSMUSG00000000001 3012.000
gene:ENSMUSG00000000003  549.500
gene:ENSMUSG00000000028 1541.881
gene:ENSMUSG00000000031 1131.651
gene:ENSMUSG00000000037 3297.491
gene:ENSMUSG00000000049  940.000
> head(txi$abundance)
                            [,1]
gene:ENSMUSG00000000001 8.148651
gene:ENSMUSG00000000003 0.000000
gene:ENSMUSG00000000028 0.331423
gene:ENSMUSG00000000031 7.391070
gene:ENSMUSG00000000037 0.259804
gene:ENSMUSG00000000049 0.007995



salmon:
Name	Length	EffectiveLength	TPM	NumReads
transcript:ENSMUST00000162897	4153	3903.000	0.044001	22.853
transcript:ENSMUST00000159265	2989	2739.000	0.000000	0.000
transcript:ENSMUST00000070533	3634	3384.000	0.064728	29.147

star:
    Geneid  Chr Start   End Strand  Length  PyMT-DTR_Replicate1_S1/mark_duplicates/PyMT-DTR_Replicate1_S1.bam
    ENSMUSG00000051951  1   3205901 3671498 -   465598  265
"""
    def get_feature_counts_prokaryotes(self, feature=None, attribute=None):
        if feature is None:
            feature = "gene"
        if attribute is None:
            attribute = "ID"

        annot = self.gff.df
        annot = annot.query("type==@feature").copy()
        names = [x[attribute] for x in annot.attributes]
        identifiers = [x[attribute] for x in annot.attributes]

        annot['identifiers'] = identifiers
        annot['names'] = names
        annot = annot.set_index("identifiers")

        ff = self.filename.split("/")[-1]
        results = f"Geneid\tChr\tStart\tEnd\tStrand\tLength\t{ff}"

        for name, length in zip(self.df.Name, self.df.Length):
            try:dd = annot.loc[name]
            except: continue
            if isinstance(dd.seqid, str):
                length2 = dd.stop - dd.start +1
                seqid = dd.seqid
                stops = dd.stop
                starts = dd.start
                strands = dd.strand
                new_name = dd['names']
            else:
                seqid = ";".join(dd.seqid.values)
                starts = ";".join([str(x) for x in dd.start.values])
                stops = ";".join([str(x) for x in dd.stop.values])
                strands = ";".join(dd.strand.values)
                length2 = (dd.stop - dd.start).sum() + len(dd.stop)
                new_name = dd['names'].values[0]

            if abs(length - length2) > 5:
                print(name, length, length2)
                raise ValueError("length in gff and quant not the same")
            NumReads = int(self.df.query("Name==@name")['NumReads'].values[0])
            if name.startswith("gene"):
                results += f"\n{name}\t{seqid}\t{starts}\t{stops}\t{strands}\t{length}\t{NumReads}"
        return results

    def save_feature_counts(self, filename, feature="gene", attribute="ID"):
        from sequana import version
        data = self.get_feature_counts(feature=feature, attribute=attribute)
        with open(filename, "w") as fout:
            fout.write(f"# Program:sequana.salmon v{version}; sequana " +
                       f"salmon -i {self.filename} -o {filename} -g {self.gff.filename}\n")
            fout.write(data)



