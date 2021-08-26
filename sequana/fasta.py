#
#  This file is part of Sequana software
#
#  Copyright (c) 2016-2021 - Sequana Development Team
#
#  Distributed under the terms of the 3-clause BSD license.
#  The full license is in the LICENSE file, distributed with this software.
#
#  website: https://github.com/sequana/sequana
#  documentation: http://sequana.readthedocs.io
#
##############################################################################
"""Utilities to manipulate FastA files"""
import os
from pysam import FastxFile
from easydev import Progress
import textwrap

from sequana.stats import N50, L50

import colorlog

logger = colorlog.getLogger(__name__)


__all__ = ["FastA"]


def is_fasta(filename):
    with open(filename, "r") as fin:
        try:
            line = fin.readline()
            assert line.startswith(">")
            line = fin.readline()
            return True
        except:  # pragma: no cover
            return False


# cannot inherit from FastxFile (no object in the API ?)
class FastA(object):
    """Class to handle FastA files


    ::

        from sequana import FastA
        f = FastA("test.fa")
        read = next(f)

        names = f.names

    """

    def __init__(self, filename, verbose=False):
        self._fasta = FastxFile(filename)
        self.filename = filename
        self._N = None

    def __iter__(self):
        return self

    def __next__(self):  # python 3
        return self.next()

    def next(self):  # python 2
        # reads 4 lines
        try:
            d = next(self._fasta)
            return d
        except KeyboardInterrupt:  # pragma: no cover
            # This should allow developers to break a loop that takes too long
            # through the reads to run forever
            self._fasta.close()
            self._fasta = FastxFile(self._fasta.filename)
        except:
            self._fasta.close()
            self._fasta = FastxFile(self._fasta.filename)
            raise StopIteration

    def __len__(self):
        if self._N is None:
            logger.info("Reading input fasta file...please wait")
            self._N = len([x for x in FastxFile(self.filename)])
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
        """Return dictionary with sequence names and lengths as keys/values"""
        return dict(zip(self.names, self.lengths))

    def format_contigs_denovo(self, output_file, len_min=500):
        """Remove contigs with sequence length below specific threshold.

        :param str output_file: output file name.
        :param int len_min: minimal length of contigs.

        Example::

            from sequana import FastA

            contigs = FastA("denovo_assembly.fasta")
            contigs.format_contigs_denovo("output.fasta", len_min=500)

        """
        # catch basename of file without extension
        project = os.path.basename(output_file).split(".")[0]
        # check if directory exist
        output_dir = os.path.dirname(output_file)
        try:
            if not os.path.exists(output_dir):  # pragma: no cover
                os.makedirs(output_dir)
        except FileNotFoundError:  # pragma: no cover
            pass

        n = 1
        with open(output_file, "w") as fp:
            for contigs in self:
                if len(contigs.sequence) < len_min:
                    break
                name = ">{}_{} {}\n".format(project, n, contigs.name)
                sequence = (
                    "\n".join(
                        [
                            contigs.sequence[i : min(i + 80, len(contigs.sequence))]
                            for i in range(0, len(contigs.sequence), 80)
                        ]
                    )
                    + "\n"
                )
                fp.write(name + sequence)
                n += 1

    def filter(self, output_filename, names_to_keep=None, names_to_exclude=None):
        """save FastA excluding or including specific sequences"""
        if names_to_exclude is None and names_to_keep is None:  # pragma: no cover
            logger.warning("No ids provided")
            return

        if names_to_exclude:
            with open(self.filename) as fin:
                with open(output_filename, "w") as fout:
                    skip = False
                    # do no use readlines. may be slower but may cause memory
                    # issue
                    for line in fin:
                        if line.startswith(">"):
                            if line[1:].split()[0] in names_to_exclude:
                                skip = True
                            else:
                                skip = False
                        if skip is False:
                            fout.write(line)
        elif names_to_keep:
            with open(self.filename) as fin:
                with open(output_filename, "w") as fout:
                    # do no use readlines. may be slower but may cause memory
                    # issue
                    skip = True
                    for line in fin:
                        if line.startswith(">"):
                            if line[1:].split()[0] in names_to_keep:
                                skip = False
                            else:
                                skip = True
                        if skip is False:
                            fout.write(line)

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
        pb = Progress(thisN)  # since we scan the entire file
        with open(output_filename, "w") as fh:
            for i, read in enumerate(fasta):
                if i in cherries:
                    fh.write(read.__str__() + "\n")
                else:
                    pass
                pb.animate(i + 1)
        return cherries

    def get_stats(self):
        """Return a dictionary with basic statistics"""
        from pylab import mean

        stats = {}
        stats["N"] = len(self.sequences)
        stats["mean_length"] = mean(self.lengths)
        stats["total_length"] = sum(self.lengths)
        stats["N50"] = N50(self.lengths)
        stats["L50"] = L50(self.lengths)
        stats["min_length"] = min(self.lengths)
        stats["max_length"] = max(self.lengths)
        return stats

    def summary(self, max_contigs=-1):
        """returns summary and print information on the stdout

        This method is used when calling sequana standalone ::

            sequana summary test.fasta
        """
        from pylab import mean, argmax

        # used by sequana summary fasta
        summary = {"number_of_contigs": len(self.sequences)}
        summary["total_contigs_length"] = sum(self.lengths)
        summary["mean_contig_length"] = mean(self.lengths)
        summary["max_contig_length"] = max(self.lengths)
        summary["min_contig_length"] = min(self.lengths)
        N = 0
        lengths = self.lengths[:]
        positions = list(range(len(lengths)))
        stats = self.get_stats()
        print("#sample_name: {}".format(self.filename))
        print("#total length: {}".format(stats["total_length"]))
        print("#N50: {}".format(stats["N50"]))
        print("#Ncontig: {}".format(stats["N"]))
        print("#L50: {}".format(stats["L50"]))
        print("#max_contig_length: {}".format(stats["max_length"]))
        print("#min_contig_length: {}".format(stats["min_length"]))
        print("#mean_contig_length: {}".format(stats["mean_length"]))

        print("contig name,length,count A,C,G,T,N")
        if max_contigs == -1:
            max_contigs = len(lengths) + 1
        while lengths and N < max_contigs:
            N += 1
            index = argmax(lengths)
            length = lengths.pop(index)
            position = positions.pop(index)
            sequence = self.sequences[position]
            name = self.names[position]
            print(
                "{},{},{},{},{},{},{}".format(
                    name,
                    length,
                    sequence.count("A"),
                    sequence.count("C"),
                    sequence.count("G"),
                    sequence.count("T"),
                    sequence.count("N"),
                )
            )

    def GC_content_sequence(self, sequence):
        """Return GC content in percentage of a sequence"""
        GC = sequence.count("G") + sequence.count("g")
        GC += sequence.count("C") + sequence.count("c")
        return GC / len(sequence) * 100

    def GC_content(self):
        """Return GC content in percentage of all sequences found in the FastA file"""
        lengths = sum(self.lengths)
        GC = 0
        for seq in self.sequences:
            GC += seq.count("G") + seq.count("g")
            GC += seq.count("C") + seq.count("c")
        return GC / lengths * 100

    def reverse_and_save(self, filename):
        """Reverse sequences and save in a file"""
        with open(filename, "w") as fout:
            for read in self:
                fout.write(
                    ">{}\t{}\n{}\n".format(read.name, read.comment, read.sequence[::-1])
                )

    def save_ctg_to_fasta(self, ctgname, outname, max_length=-1):
        """Select a contig and save in a file"""
        index = self.names.index(ctgname)
        with open("{}.fa".format(outname), "w") as fout:

            if max_length == -1:
                fout.write(">{}\n{}".format(outname, self.sequences[index]))
            else:
                fout.write(
                    ">{}\n{}".format(outname, self.sequences[index][0:max_length])
                )

    def to_fasta(self, outfile, width=80):
        """Save the input FastA file into a new file

        The interest of this method is to wrap the sequence into 80 characters.
        This is useful if the input file is not formatted correctly.

        """
        with open(outfile, "w") as fout:
            for name, comment, seq in zip(self.names, self.comments, self.sequences):
                seq = "\n".join(textwrap.wrap(seq, width))
                if comment is None:
                    fout.write(f">{name}\n{seq}\n")
                else:
                    fout.write(f">{name}\t{comment}\n{seq}\n")

    def to_igv_chrom_size(self, output):
        """Create a IGV file storing chromosomes and their sizes"""
        data = self.get_lengths_as_dict()
        with open(output, "w") as fout:
            for k, v in data.items():
                fout.write("{}\t{}\n".format(k, v))

    def save_collapsed_fasta(self, outfile, ctgname, width=80, comment=None):
        """Concatenate all contigs and save results"""
        with open(outfile, "w") as fout:
            data = "".join(self.sequences)
            seq = "\n".join(textwrap.wrap(data, width))
            if comment is None:
                fout.write(f">{ctgname}\n{seq}\n")
            else:
                fout.write(f">{ctgname}\t{comment}\n{seq}\n")
