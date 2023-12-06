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
"""Utilities to manipulate and find codons"""


class Codon:
    """Utilities to manipulate codons

    The codon contains hard-coded set of start and stop codons (bacteria) for
    strand plus and minus. Adapt to your needs for other organisms. Based on the
    scan of Methanothermobacter thermautotrophicus bacteria.

    ::

        from sequana import Codon
        c = Codon()
        c.start_codons['+']

    """

    codons = {
        "start": {
            "+": frozenset({"ATG", "TTG", "GTG"}),
            "-": frozenset({"CAT", "CAA", "CAC"}),
        },
        "stop": {
            "+": frozenset({"TAG", "TGA", "TAA"}),
            "-": frozenset({"TTA", "TCA", "CTA"}),
        },
    }

    def __init__(self):
        pass

    def get_codons_from_fasta_and_gff(self, fasta, gff):
        raise NotImplementedError

    def find_start_codon_position(self, sequence, position, strand, max_shift=10000):
        """Return starting position and string of closest start codon to a given position

        **The starting position is on the 5-3 prime direction (see later)**

        :param str sequence:
        :param int position: 0-base position
        :param str strand: '+' or '-'

        The search starts at the given position, then +1 base, then -1 base,
        then +2, -2, +3, etc

        Here, we start at position 2 (letter G), then shift by +1 and find the ATG
        string.
        ::

            >>> from sequana import Codon
            >>> c = Codon()
            >>> c.find_start_codon_position("ATGATGATG", 2, "+")
            (3, 'ATG')

        whereas starting at position 1, a shift or +1 (position 2 ) does not hit a start codon.
        Next, a shift of -1 (position 0) hits the ATG start codon.::

            >>> c.find_start_codon_position("ATGATGATG", 1, "+")
            (3, 'ATG')

        On strand -, the start codon goes from right to left. So, in the following example,
        the CAT start codon (reverse complement of ATG) is found at position 3. Developers must take
        into account a +3 shift if needed::

            >>> c.find_start_codon_position("AAACAT", 3, "-")
            (3, 'CAT')

            >>> c.find_start_codon_position("AAACATCAT", 8, "-")

        """
        assert strand in ("+", "-")
        codons = self.codons["start"][strand]
        return self._search_codons(sequence, position, strand, max_shift, codons)

    def find_stop_codon_position(self, sequence, position, strand, max_shift=10000):
        """Return position and string of closest stop codon to a given position

        :param str sequence:
        :param int position: 0-base position
        :param str strand: + or -

        See :meth:`find_start_codon_position` for details.

        Only difference is that the search is based on stop codons rather than start codons.

        ::

            >>> from sequana import Codon
            >>> c = Codon()
            >>> c.find_stop_codon_position("ATGACCCC", 2, "+")
            (1, 'TGA')

        """
        assert strand in ("+", "-")
        codons = self.codons["stop"][strand]
        return self._search_codons(sequence, position, strand, max_shift, codons)

    def _search_codons(self, sequence, position, strand, max_shift, codons):
        max_shift = max(max_shift, len(sequence))

        # We alternate the starting position on each side starting with a shift of 0,
        # then 1, then -1, 2, -2 and so on
        found = False
        shift = 0

        # on strand +, we start the shift on the right (inside the gene). On strand -, we
        # shift to the left first (0,-1,1,-2,2).
        side = 1 if strand == "+" else -1

        while not found:
            # if we reach the position 0, no need to search for the codon
            ps = position + shift
            if ps < 0:
                pass
            elif sequence[ps : ps + 3] in codons:
                return ps, sequence[ps : ps + 3]

            if side == 1:
                shift = -shift + 1
                if shift > max_shift:  # pragma: no cove #pragma: no cover
                    break
            else:
                shift *= -1
            side *= -1
