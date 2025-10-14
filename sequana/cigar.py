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
import re
from collections import defaultdict

"""
Note: could use pysam most probably to improve the speed.
"""
import colorlog

logger = colorlog.getLogger(__name__)


class Cigar:
    """A class to handle CIGAR strings from BAM files.

    .. doctest::

        >>> from sequana.cigar import Cigar
        >>> c = Cigar("2S30M1I")
        >>> len(c)
        33

        >>> c = Cigar("1S1S1S1S")
        >>> c.compress()
        >>> c.cigarstring
        '4S'


    Possible CIGAR types are:

    - "M" : alignment match
    - "I" : insertion to the reference
    - "D" : deletion from the reference
    - "N" : skipped region from the reference
    - "S" : soft clipping (clipped sequence present in seq)
    - "H" : hard clipping (sequence NOT present)
    - "P" : padding (silent deletion from padded reference)
    - "=" : sequence match
    - "X" : sequence mismatched
    - "B" : back (rare) (could be also NM ?)

    !!! BWA MEM get_cigar_stats returns list with 11 items
    Last item is
    !!! what is the difference between M and = ???
    Last item is I + S + X
    !!! dans BWA, mismatch (X) not provided... should be deduced from last item - I - S

    .. note:: the length of the query sequence based on the CIGAR is calculated
        by adding the M, I, S, =, or X and other operations are ignored.
        source: https://stackoverflow.com/questions/39710796/infer-the-length-of-a-sequence-using-the-cigar/39812985#39812985

    :reference: https://github.com/samtools/htslib/blob/develop/htslib/sam.h
    """

    __slots__ = ("_cigarstring",)
    pattern = re.compile(r"(\d+)([A-Za-z])?")
    types = "MIDNSHP=XB"

    def __init__(self, cigarstring: str):
        """.. rubric:: Constructor

        :param str cigarstring: the CIGAR string.

        .. note:: the input CIGAR string validity is not checked.
            If an unknown type is found, it is ignored generally.
            For instance, the length of 1S100Y is 1 since Y is not correct.

        """
        if not isinstance(cigarstring, str):
            raise TypeError("Cigar string must be a string.")
        #: the CIGAR string attribute
        self._cigarstring = cigarstring

    @property
    def cigarstring(self):
        return self._cigarstring

    def __str__(self):
        return self._cigarstring

    def __repr__(self):
        return "Cigar( {} )".format(self._cigarstring)

    def __len__(self):
        return sum([y for x, y in self._decompose() if x in "MIS=X"])

    def _decompose(self):
        for num, op in self.pattern.findall(self._cigarstring):
            if op and op in self.types:
                yield op, int(num)

    def as_sequence(self):
        return "".join(op * count for op, count in self._decompose())

    def as_dict(self):
        """Return cigar types and their count

        :return: dictionary

        Note that repeated types are added::

            >>> c = Cigar('1S2M1S')
            >>> c.as_dict()
            {"S":2,"M":2}

        """
        # !! here, we have to make sure that  duplicated letters are summed up

        d = defaultdict(int)
        for op, count in self._decompose():
            d[op] += count
        return d

    def as_tuple(self):
        """Decompose the cigar string into tuples keeping track of repeated types

        :return: tuple

        .. doctest::

            >>> from sequana import Cigar
            >>> c = Cigar("1S2M1S")
            >>> c.as_tuple()
            (('S', 1), ('M', 2), ('S', 1))

        """
        return tuple(self._decompose())

    def compress(self):
        """1S1S should become 2S. inplace modification"""
        data = list(self._decompose())
        if not data:
            return

        compressed = [data[0]]
        for op, count in data[1:]:
            last_op, last_count = compressed[-1]
            if op == last_op:
                compressed[-1] = (op, last_count + count)
            else:
                compressed.append((op, count))
        self._cigarstring = "".join(f"{count}{op}" for op, count in compressed)

    def stats(self):
        """Returns number of occurence for each type found in :attr:`types`

        ::

            >>> c = Cigar("1S2M1S")
            >>> c.stats()
            [2, 0, 0, 0, 2, 0, 0, 0, 0, 0]

        """
        counts = [0] * len(self.types)
        d = self.as_dict()
        for op, val in d.items():
            idx = self.types.find(op)
            if idx != -1:
                counts[idx] = val
        return counts

    def get_query_length(self) -> int:
        """Return length consumed by the query (read)."""
        return sum(count for op, count in self._decompose() if op in "MIS=X")

    def get_reference_length(self) -> int:
        """Return length consumed by the reference."""
        return sum(count for op, count in self._decompose() if op in "MDN=X")


def fetch_exon(chrom, start, cigar):
    chrom_start = start
    exon_bound = []
    for c, size in cigar:
        if c == 0:
            exon_bound.append((chrom, chrom_start, chrom_start + size))
            chrom_start += size
        elif c == 1:
            continue
        elif c == 2:
            chrom_start += size
        elif c == 3:
            chrom_start += size
        elif c == 4:  # FIXME do we want to include this in the exon
            chrom_start += size
        else:
            continue
    return exon_bound


def fetch_intron(chrom, start, cigar):
    # equivalence:
    # c = 0 -> M
    # c = 1 -> I
    # c = 2 -> D
    # c = 3 -> N gap/intron
    # c = 4 -> S soft clipping
    chrom_start = start
    intron_bound = []
    for c, size in cigar:
        if c == 0:
            chrom_start += size
        elif c == 1:
            continue
        elif c == 2:
            chrom_start += size
        elif c == 3:
            intron_bound.append((chrom, chrom_start, chrom_start + size))
            chrom_start += size
        elif c == 4:  # not including soft clipping in the intron
            continue
        else:
            continue
    return intron_bound


def fetch_clip(chrom, start, cigar):
    chrom_start = start
    clip_bound = []
    for c, size in cigar:
        if c == 0:
            chrom_start += size
        elif c == 1:
            continue
        elif c == 2:
            chrom_start += size
        elif c == 3:
            chrom_start += size
        elif c == 4:
            clip_bound.append((chrom, chrom_start, chrom_start + size))
            chrom_start += size
        else:
            continue
    return clip_bound


def fetch_deletion(chrom, start, cigar):
    chrom_start = start
    deletion_bound = []
    for c, size in cigar:
        if c == 0:
            chrom_start += size
        elif c == 1:
            continue
        elif c == 2:
            deletion_bound.append((chrom, chrom_start, chrom_start + size))
            chrom_start += size
        elif c == 3:
            chrom_start += size
        elif c == 4:
            chrom_start += size
        else:
            continue
    return deletion_bound


def fetch_insertion(chrom, start, cigar):
    # NOTE that the returned insertions are stored as
    # chrom, start, size rather than chrom, start, end in other fetchers
    # functions
    chrom_start = start
    insertion_bound = []
    for c, size in cigar:
        if c == 0:
            chrom_start += size
        elif c == 1:
            # See note above
            insertion_bound.append((chrom, chrom_start, size))
            continue
        elif c == 2:
            chrom_start += size
        elif c == 3:
            chrom_start += size
        elif c == 4:
            chrom_start += size
        else:
            continue
    return insertion_bound
