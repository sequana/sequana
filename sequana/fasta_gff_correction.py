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
"""Utilities to manipulate FastA files"""
import textwrap
from collections import defaultdict

import colorlog
from pysam import FastxFile

from sequana import GFF3, Codon, FastA, VCF_freebayes
from sequana.lazy import pandas as pd

logger = colorlog.getLogger(__name__)


__all__ = ["FastaGFFCorrection"]


class FastaGFFCorrection:
    """Class to correct FastA reference and its GFF given a set of variants

    ::
        f = FastaGFFCorrection("input.fa", "input.vcf")
        f.fix_and_save_fasta("new.fa")
        f.fix_and_save_gff("new.fa", "input.gff3", "new.gff3")

    You must call :meth:`save_fasta` before :meth:`save_gff` to compute
    the shifts of positions found in the VCF file.

    The GFF corrections are made on the start and stop positions of all features
    found in the GFF. In addition, gene and CDS start and stop positions are also
    corrected for the change of start/stop codon implied by the variants.


    """

    def __init__(self, fasta_file, vcf_file):
        self._fasta_file = fasta_file

        self.vcf = VCF_freebayes(vcf_file)

        self.shifts = None
        self.positions = None

        self.codon = Codon()

        # parameter used in find_start_codon_position method.
        self._start_codon_max_shift = 10000

    def fix_and_save_fasta(self, outfile):
        """

        We consider that the VCF file contains variants that have already been filtered.

        .. todo:: handle multiple-variants at a given position
        """

        # initiate the fasta file
        self._fasta = FastxFile(self._fasta_file)

        # read the variants once for all
        variants = self.vcf.get_variants()

        # and build a convenient dataframe
        df_vcf = pd.DataFrame([v.resume for v in variants])

        # let us make sure the positions are integers
        df_vcf["position"] = [int(x) for x in df_vcf.position]

        # what chrom names do we have ?
        chrom_names = df_vcf["chr"].unique()

        # let us open the output FastA file handler
        with open(outfile, "w") as fout:
            # we will keep track of shifts due to INDELs and their positions
            # for each chromosome. This information will be required if we also
            # want to correct the GFF file.
            shifts = {k: [] for k in chrom_names}
            positions = {k: [] for k in chrom_names}

            # let us scan each input sequence in the input FastA file
            for entry in self._fasta:
                # retrieve name/comment/sequence
                sequence, chrname, comment = entry.sequence, entry.name, entry.comment

                # initiate the new sequence with variants
                new_sequence = ""

                # we will use a cursor to scan the original sequence
                cursor = 0

                diff = 0
                for pos, ref, alt in df_vcf.query("chr == @chrname")[["position", "reference", "alternative"]].values:
                    # accumulate bases before the variant and add the alternate
                    new_sequence += sequence[cursor : pos - 1] + alt

                    # we update the cursor to scan the original sequence
                    cursor = pos - 1 + len(ref)

                    assert sequence[pos - 1 : pos - 1 + len(ref)] == ref

                    if len(ref) != len(alt):
                        diff += len(alt) - len(ref)
                        shifts[chrname].append(diff)
                        positions[chrname].append(pos)

                # do not forget the final slice after last variant
                new_sequence += sequence[cursor:]

                fout.write(f">{chrname}\t{comment}\n")
                for line in textwrap.wrap(new_sequence, width=70):
                    fout.write(line + "\n")

        self.shifts = shifts
        self.positions = positions
        return {"shifts": shifts, "positions": positions}

    def _get_shift(self, position, chrom_name):
        """Utility method to get the shift created by INDELs at a given position


        This functions returns the sum of Insertion minus Deletions on the LHS of
        a given position
        """
        shifts = self.shifts[chrom_name]
        positions = self.positions[chrom_name]

        if position < positions[0]:
            return 0
        else:
            # identify the closest INDELs just before
            # got forward and stop once we are before the position
            for pos, shift in zip(positions[::-1], shifts[::-1]):
                if pos <= position:
                    return shift

    def fix_and_save_gff(self, gff_infile, fasta_corrected, gff_outfile, maxshift=8000):
        """Corrects the start/stop codons of input GFF given corrected fasta

        You must call :meth:`fix_and_save_fasta` first.

        .. todo:: multiple contig is not yet implemented"""
        # initiate the original fasta file
        fasta = FastxFile(self._fasta_file)
        contig = next(fasta)
        L = len(contig.sequence)

        # reads the corrected fasta
        corfasta = FastxFile(fasta_corrected)
        contig = next(corfasta)
        chrom_name = contig.name
        corseq = contig.sequence
        newL = len(corseq)

        with open(gff_infile, "r") as fin, open(gff_outfile, "w") as fout:
            for line in fin:
                # we can skip empty lines
                if not line.strip():
                    continue

                if line.startswith("#"):
                    # first commented lines usually contains the sequence
                    # length. may also contain the sequence name
                    line = line.replace(str(L), str(newL))
                    fout.write(f"{line}")
                else:
                    # process line, splitting items
                    items = line.split("\t")

                    # regions just need to update the length of the region
                    if items[2] == "region":
                        items[4] = str(newL)
                        line = "\t".join(items)
                        fout.write(f"{line}")
                    else:
                        # we update the starting and ending positions of all other features
                        start, end = int(items[3]), int(items[4])

                        S1 = self._get_shift(start, chrom_name)
                        S2 = self._get_shift(end, chrom_name)

                        start += S1
                        end += S2

                        # We could stop here but gene/CDS need to be check carefully since
                        # their starting/ending codon may not be correct anymore due
                        # to the shift (not a modulo 3 anymore), or variant that may happen
                        # exactly at the start/stop codon position. In principle the start codon
                        # should be correct (except if a variant occured) but stop codon will differ
                        # if an INDEL exists in between
                        if items[2] in ["gene", "CDS"]:
                            strand = items[6]
                            if strand == "+":
                                # given a gene/CDS, first retrieve the start codon position
                                # on the corrected sequence at start - 1 (-1 because python
                                # uses 0-base convetion.
                                start = self.find_start_codon_position(corseq, start - 1, strand="+") + 1

                                assert corseq[start - 1 : start - 1 + 3] in self.codon.codons["start"]["+"]

                                # We now scan the CDS/gene from start to the end searching for stop codon
                                # until we reach the end position.
                                stop_codon_position = 0
                                for cursor in range(start, end + 1, 3):
                                    codon = corseq[cursor - 1 : cursor - 1 + 3]
                                    if codon in self.codon.codons["stop"]["+"]:
                                        stop_codon_position = cursor

                                # we now have the closest stop codon to the end position
                                # we keep the delta with respect to the current known end position
                                delta = end - stop_codon_position

                                # and see if we can get a closer one after the end position
                                for cursor in range(end + 1, end + maxshift, 3):
                                    codon = corseq[cursor - 1 : cursor - 1 + 3]
                                    if codon in self.codon.codons["stop"]["+"]:
                                        stop_codon_position = cursor
                                        break

                                    if cursor - end > delta and stop_codon_position != 0:
                                        break

                                # we need to check whether we found something
                                if stop_codon_position == 0:
                                    raise Exception("stop_codon_position no found")
                                items[3] = str(start)
                                items[4] = str(stop_codon_position + 2)

                            elif strand == "-":
                                # given a gene/CDS, first retrieve the start codon position
                                # on the corrected sequence at start - 1 (-1 because python
                                # uses 0-base convetion.
                                end = self.find_start_codon_position(corseq, end - 3 - 1, strand="-")
                                end += 3 + 1  # 3 for the end of the codon +1 for 1-base convention

                                assert corseq[end - 3 - 1 : end - 1] in self.codon.codons["start"]["-"]

                                # We now scan the CDS/gene from start to the end searching for stop codon
                                # until we reach the end position.
                                stop_codon_position = 0
                                for cursor in range(end, start - 1, -3):
                                    codon = corseq[cursor - 3 - 1 : cursor - 1]
                                    if codon in self.codon.codons["stop"]["-"]:
                                        stop_codon_position = cursor - 3

                                # we now have the closest stop codon to the end position
                                # we keep the delta with respect to the current known end position
                                delta = stop_codon_position - start - 1

                                # and see if we can get a closer one after the end position
                                for cursor in range(start, start - maxshift - 1, -3):
                                    codon = corseq[cursor - 3 - 1 : cursor - 1]
                                    if codon in self.codon.codons["stop"]["-"]:
                                        stop_codon_position = cursor - 3
                                        break

                                    if start - 1 - cursor > delta and stop_codon_position != 0:
                                        break

                                # we need to check whether we found something
                                if stop_codon_position == 0:
                                    raise Exception("stop_codon_position not found")
                                items[3] = str(stop_codon_position)
                                items[4] = str(end - 1)
                            else:
                                raise ValueError(f"strand must be + or -. Found {strand}")
                        line = "\t".join(items)
                        fout.write(f"{line}")

    def find_start_codon_position(self, sequence, position, strand):
        """return position of start codon"""

        try:
            position, codon = self.codon.find_start_codon_position(
                sequence, position, strand, max_shift=self._start_codon_max_shift
            )
            return position
        except TypeError:
            logger.warning(f"Could not find start codon at position {position}. keep given position")
            return position

    def get_all_start_codons(self, fasta, gff, chr_name, strand, feature="CDS"):
        """

        On strand +, start codon are ATG, TTG, GTG
        On strand -, start codon (from 3-5 prime) are CAT, CAA, CAC (reverse complement
        of start codon on strand+).


        """

        gff = GFF3(gff)
        fasta = FastA(fasta)

        for ctg in fasta:
            if ctg.name == chr_name:
                seq = ctg.sequence
                break

        d = defaultdict(int)

        if strand == "+":
            for start in gff.df.query("seqid==@ctg.name and genetic_type==@feature and strand=='+'").start:
                d[ctg.sequence[start - 1 : start + 3 - 1]] += 1
        else:
            # for start in ...stop is not a mistake
            for start in gff.df.query("seqid==@ctg.name and genetic_type==@feature and strand=='-'").stop:
                d[ctg.sequence[start - 3 : start]] += 1
        return d

    def get_all_stop_codons(self, fasta, gff, chr_name, strand, feature="CDS"):
        """

        On strand +, stop codons are TAG, TGA, TAA
        On strand -, stop codons (from 3-5 prime) are TTA, TCA, CTA (reverse complement
        of start codon on strand+).

        """

        gff = GFF3(gff)
        fasta = FastA(fasta)

        for ctg in fasta:
            if ctg.name == chr_name:
                seq = ctg.sequence
                break

        d = defaultdict(int)

        if strand == "+":
            for stop in gff.df.query("seqid==@ctg.name and genetic_type==@feature and strand=='+'").stop:
                d[ctg.sequence[stop - 3 : stop]] += 1
        else:
            # for stop in ...start is not a mistake
            for stop in gff.df.query("seqid==@ctg.name and genetic_type==@feature and strand=='-'").start:
                d[ctg.sequence[stop - 1 : stop - 1 + 3]] += 1
        return d
