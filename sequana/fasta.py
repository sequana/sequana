# -*- coding: utf-8 -*-
#
#  This file is part of Sequana software
#
#  Copyright (c) 2016 - Sequana Development Team
#
#  File author(s):
#      Thomas Cokelaer <thomas.cokelaer@pasteur.fr>
#      Dimitri Desvillechabrol <dimitri.desvillechabrol@pasteur.fr>, 
#          <d.desvillechabrol@gmail.com>
#
#  Distributed under the terms of the 3-clause BSD license.
#  The full license is in the LICENSE file, distributed with this software.
#
#  website: https://github.com/sequana/sequana
#  documentation: http://sequana.readthedocs.io
#
##############################################################################
"""Utilities to manipulate FASTQ and Reads"""
import os
from pysam import FastxFile
from easydev import Progress

from sequana import logger
logger.name = __name__


__all__ = ["FastA"]


def is_fasta(filename):
    with open(filename, "r") as fin:
        try:
            line = fin.readline()
            assert line.startswith(">")
            line = fin.readline()
            return True
        except:
            return False


# cannot inherit from FastxFile (no object in the API ?)
class FastA(object):
    """Class to handle FastA files. Cannot be compressed


    """
    def __init__(self, filename, verbose=False):
        if filename.endswith(".gz"):
            raise ValueError("Must be decompressed.")
        self._fasta = FastxFile(filename)
        self.filename = filename
        logger.info("Reading input fasta file...please wait") 
        self._N = len([x for x in FastxFile(filename)])

    def __iter__(self):
        return self

    def __next__(self): # python 3
        return self.next()

    def next(self): # python 2
        # reads 4 lines
        try:
            d = next(self._fasta)
            return d
        except KeyboardInterrupt:
            # This should allow developers to break a loop that takes too long
            # through the reads to run forever
            self._fasta.close()
            self._fasta = FastxFile(self._fasta.filename)
        except:
            self._fasta.close()
            self._fasta = FastxFile(self._fasta.filename)
            raise StopIteration
        return d

    def __len__(self):
        return self._N

    def _get_names(self):
        return [this.name for this in self]
    names = property(_get_names)

    def _get_sequences(self):
        return [this.sequence for this in self]
    sequences = property(_get_sequences)

    def _get_comment(self):
        return [this.comment for this in self]
    comments = property(_get_comment)

    def _get_lengths(self):
        return [len(this.sequence) for this in self]
    lengths = property(_get_lengths)

    def get_lengths_as_dict(self):
        return dict(zip(self.names, self.lengths))

    def format_contigs_denovo(self, output_file, len_min=500):
        """Replace NODE with the project name and remove contigs with a length 
        lower than len_min.

        :param str output_file: output file name.
        :param int len_min: minimal length of contigs.

        Example:

            from sequana import FastA

            contigs = FastA("denovo_assembly.fasta")
            contigs.format_contigs_denovo("path/to/file.fasta", len_min=500)

        Results are stored in "path/to/file.fasta".
        """
        # catch basename of file without extension
        project = os.path.basename(output_file).split(".")[0]
        # check if directory exist
        output_dir = os.path.dirname(output_file)
        try:
            if not os.path.exists(output_dir):
                os.makedirs(output_dir)
        except FileNotFoundError:
            pass

        n = 1
        with open(output_file, "w") as fp:
            for contigs in self:
                if len(contigs.sequence) < len_min:
                    break
                name = ">{}_{} {}\n".format(project, n, contigs.name)
                sequence = "\n".join([contigs.sequence[i:min(i+80, 
                    len(contigs.sequence))] for i in range(0, 
                    len(contigs.sequence), 80)]) + "\n"
                fp.write(name + sequence)
                n += 1

    def select_random_reads(self, N=None, output_filename="random.fasta"):
        """Select random reads and save in a file

        :param int N: number of random unique reads to select
            should provide a number but a list can be used as well.
        :param str output_filename:
        """
        import numpy as np
        thisN = len(self)
        if isinstance(N, int):
            if N > thisN:
                N = thisN
            # create random set of reads to pick up
            cherries = list(range(thisN))
            np.random.shuffle(cherries)
            # cast to set for efficient iteration
            cherries = set(cherries[0:N])
        elif isinstance(N, set):
            cherries = N
        elif isinstance(N, list):
            cherries = set(N)
        fasta = FastxFile(self.filename)
        pb = Progress(thisN) # since we scan the entire file
        with open(output_filename, "w") as fh:
            for i, read in enumerate(fasta):
                if i in cherries:
                    fh.write(read.__str__() + "\n")
                else:
                    pass
                pb.animate(i+1)
        return cherries

    def get_stats(self):
        from pylab import mean
        stats = {}
        stats["N"] = len(self.sequences)
        stats["mean_length"] = mean(self.lengths)
        return stats

    def reverse_and_save(self, filename):
        with open(filename, "w") as fout:
            for read in self:
                fout.write(">{}\t{}\n{}\n".format(read.name, read.comment,
                    read.sequence[::-1]))

    def save_ctg_to_fasta(self, ctgname, outname):
        index = self.names.index(ctgname)
        with open("{}.fa".format(outname), "w") as fout:
            fout.write(">{}\n{}".format(outname, self.sequences[index]))

    def to_fasta(self, outfile, width=80):
        """Save the input FastA file into a new file

        The interest of this method is to wrap the sequence into 80 characters.
        This is useful if the input file is not formatted correctly.

        """
        with open(outfile, "w") as fout:
            for name,comment,seq in zip(self.names, self.comments, self.sequences):
                import textwrap
                seq = "\n".join(textwrap.wrap(seq, width))
                if comment is None:
                    fout.write(">{}\n{}\n".format(name, seq))
                else:
                    fout.write(">{}\t{}\n{}\n".format(name, comment, seq))

    def to_igv_chrom_size(self, output):
        data = self.get_lengths_as_dict()
        with open(output, "w") as fout:
        
            for k,v in data.items():
                fout.write("{}\t{}\n".format(k, v))











