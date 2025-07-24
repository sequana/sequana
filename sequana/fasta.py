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
import textwrap
from collections import defaultdict

import colorlog
import tqdm

from sequana.lazy import numpy as np
from sequana.lazy import pylab, pysam
from sequana.stats import L50, N50

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
class FastA:
    """Class to handle FastA files

    Works like an iterator::

        from sequana import FastA
        f = FastA("test.fa")
        read = next(f)

    names and sequences can be accessed with attributes::

        f.names
        f.sequences


    """

    def __init__(self, filename, verbose=False):
        self._fasta = pysam.FastaFile(filename)
        self._fastx = pysam.FastxFile(filename)
        self.filename = filename
        self._N = None

    def __iter__(self):
        return self

    def __next__(self):  # python 3
        return self.next()

    def next(self):  # python 2
        try:
            d = next(self._fastx)
            return d
        except KeyboardInterrupt:  # pragma: no cover
            # This should allow developers to break a loop that takes too long
            # through the reads to run forever
            self._fastx.close()
            self._fastx = pysam.FastxFile(self.filename)
        except:
            self._fastx.close()
            self._fastx = pysam.FastxFile(self.filename)
            raise StopIteration

    def __getitem__(self, name):
        return self._fasta[name]

    def __len__(self):
        return len(self._fasta)

    def _get_names(self):
        return self._fasta.references

    names = property(_get_names)

    def _get_sequences(self):
        return [self._fasta.fetch(name) for name in self.names]

    sequences = property(_get_sequences)

    def _get_comment(self):
        self._fastx = pysam.FastxFile(self.filename)
        return [this.comment for this in self._fastx]

    comments = property(_get_comment)

    def get_cumulative_sum(self, mode="mixed", exclude=[]):
        """Compute the cumulative sum of values from a dictionary, sorted by name.

        This method returns two lists:
        - A list of names sorted according to the specified mode.
        - A list of cumulative sums corresponding to the lengths of the sorted names.

        Sorting behavior:
        - If `mode="mixed"` (default), names containing numbers are sorted naturally,
          meaning numerical values are sorted as integers while preserving non-numeric strings.
        - If `mode="alphanum"`, all names are treated as strings and sorted lexicographically.

        Parameters:
        - mode (str): Sorting mode, either `"mixed"` (natural sorting) or `"alphanum"` (lexicographic).
        - exclude (list): List of names to exclude from processing.

        Returns:
            - tuple: (sorted_names, cumulative_sums)
                - sorted_names (list): Names sorted based on the selected mode.
                - cumulative_sums (list): Cumulative sum of corresponding values.

        Example:

        If input is `['1', 'maxi', '10', '2']`, mixed mode returns `['1', '2', '10', 'maxi']`,
        ensuring `['1', 10, 2, 'maxi']` is correctly ordered as `[1, 2, '10', 'maxi']`.
        """

        assert mode in ["mixed", "alphanum"]

        if mode == "mixed":
            sorted_names = self.sorted_mixed_names
        else:
            sorted_names = self.sorted_names

        from itertools import accumulate

        lengths = self.get_lengths_as_dict()
        return sorted_names, list(accumulate([lengths[name] for name in sorted_names]))

    def _get_sorted_mixed_names(self):
        return sorted(self.names, key=lambda x: (isinstance(x, str) and not x.isdigit(), int(x) if x.isdigit() else x))

    sorted_mixed_names = property(_get_sorted_mixed_names)

    def _get_sorted_names(self):
        return sorted(self.names, key=lambda x: (isinstance(x, str) and not x.isdigit(), int(x) if x.isdigit() else x))

    sorted_names = property(_get_sorted_names)

    def _get_lengths(self):
        return self._fasta.lengths

    lengths = property(_get_lengths)

    def get_lengths_as_dict(self):
        """Return dictionary with sequence names and lengths as keys/values"""
        return dict(zip(self.names, self.lengths))

    def explode(self, outdir="."):
        """extract sequences from one fasta file and save them into individual files"""
        with open(self.filename, "r") as fin:
            for line in fin.readlines():
                if line.startswith(">"):
                    # ignore the comment and > character and use it as the
                    # filename
                    name = line.split()[0][1:]
                    try:
                        # if a file was already open, let us close it
                        fout.close()
                    except NameError:
                        pass
                    finally:
                        fout = open(f"{outdir}/{name}.fasta", "w")
                    fout.write(line)
                else:
                    fout.write(line)
        # need to close the last file
        fout.close()

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

        comments = self.comments
        with open(output_filename, "w") as fh:
            for i in cherries:
                name = self.names[i]
                seq = self._fasta.fetch(self.names[i])
                comment = comments[i]
                fh.write(f">{name}\t{comment}\n{seq}\n")
        return cherries

    def get_stats(self):
        """Return a dictionary with basic statistics

        N the number of contigs, the N50 and L50, the min/max/mean
        contig lengths and total length.
        """
        stats = {}
        stats["N"] = len(self.sequences)
        stats["mean_length"] = pylab.mean(self.lengths)
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

        # used by sequana summary fasta
        summary = {"number_of_contigs": len(self.sequences)}
        summary["total_contigs_length"] = sum(self.lengths)
        summary["mean_contig_length"] = pylab.mean(self.lengths)
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
            index = pylab.argmax(lengths)
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
        for name in self.names:
            seq = self._fasta.fetch(name)
            GC += seq.count("G") + seq.count("g")
            GC += seq.count("C") + seq.count("c")
        return GC / lengths * 100

    def reverse_and_save(self, filename):
        """Reverse sequences and save in a file"""
        with open(filename, "w") as fout:
            for read in self:
                fout.write(">{}\t{}\n{}\n".format(read.name, read.comment, read.sequence[::-1]))

    def save_ctg_to_fasta(self, ctgname, outname, max_length=-1):
        """Select a contig and save in a file"""
        with open("{}.fa".format(outname), "w") as fout:
            if max_length == -1:
                fout.write(">{}\n{}".format(outname, self._fasta.fetch(ctgname)))
            else:
                fout.write(">{}\n{}".format(outname, self._fasta.fetch(ctgname)[0:max_length]))

    def to_fasta(self, outfile, width=80, sorting="natsort"):
        """Save the input FastA file into a new file

        The interest of this method is to wrap the sequence into 80 characters.
        This is useful if the input file is not formatted correctly.

        """

        if sorting == "natsort":
            import natsort

            names = natsort.natsorted(self.names)
        else:
            names = self.names

        with open(outfile, "w") as fout:
            for name in names:
                comment = self.comments[self.names.index(name)]
                # fetch sequence and wrap it
                seq = self._fasta.fetch(name)
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

    def find_gaps(self):
        """Identify NNNNs in data

        returns a dictionary. keys are the chromosomes' names
        values is a list. the first item is the number of Ns. the next items are the gaps' positions
        """
        results = defaultdict(list)
        for i, seq in enumerate(self.sequences):
            count = 0
            positions = [0]
            for pos, x in enumerate(seq):
                if x == "N":
                    count += 1
                    positions.append(pos)
            if count:
                name = self.names[i]
                results[name].append(count)
                for i, pos in enumerate(positions[1:]):
                    if positions[i] - pos == -1:
                        pass
                    else:
                        results[name].append(pos)
        return results

    def print_sequence_region_for_gff(self):
        for name in sorted(self.names):
            L = len(self.sequences[self.names.index(name)])
            # the 2 # are important e.g. for snpeff
            print(f"##sequence-region\t{name}\t{1}\t{L}")
