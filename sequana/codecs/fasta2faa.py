#
#  This file is part of Sequana software
#
#  Copyright (c) 2026 - Sequana Development Team
#
#  Distributed under the terms of the 3-clause BSD license.
#  The full license is in the LICENSE file, distributed with this software.
#
#  website: https://github.com/sequana/sequana
#  documentation: http://sequana.readthedocs.io
#
##############################################################################
"""Convert a nucleotide FASTA file into a protein FASTA (FAA) file.

Naive, frame-0, single-strand translation of each sequence (no ORF
detection, no reverse-complement, no alternative reading frames): this is
the same behaviour as `bioconvert`'s ``fasta2faa`` converter, reimplemented
here (using sequana's own pysam-backed :class:`~sequana.fasta.FastA` reader)
because this conversion is used often enough in Sequana workflows to warrant
its own dependency-free entry point rather than requiring bioconvert.

Stop codons (``TAA``, ``TAG``, ``TGA``) are translated to ``_`` (matching
bioconvert's convention); any other codon not in the standard table (partial
trailing codon, ambiguous/non-ACGT bases) is translated to ``X``.
"""
import textwrap

from sequana.fasta import FastA

__all__ = ["fasta2faa", "translate_codon_naive", "CODON_TABLE"]

# Standard genetic code, DNA codon -> single-letter amino acid.
# Stop codons map to "_" (matches bioconvert's fasta2faa convention).
CODON_TABLE = {
    "ATA": "I", "ATC": "I", "ATT": "I", "ATG": "M",
    "ACA": "T", "ACC": "T", "ACG": "T", "ACT": "T",
    "AAC": "N", "AAT": "N", "AAA": "K", "AAG": "K",
    "AGC": "S", "AGT": "S", "AGA": "R", "AGG": "R",
    "CTA": "L", "CTC": "L", "CTG": "L", "CTT": "L",
    "CCA": "P", "CCC": "P", "CCG": "P", "CCT": "P",
    "CAC": "H", "CAT": "H", "CAA": "Q", "CAG": "Q",
    "CGA": "R", "CGC": "R", "CGG": "R", "CGT": "R",
    "GTA": "V", "GTC": "V", "GTG": "V", "GTT": "V",
    "GCA": "A", "GCC": "A", "GCG": "A", "GCT": "A",
    "GAC": "D", "GAT": "D", "GAA": "E", "GAG": "E",
    "GGA": "G", "GGC": "G", "GGG": "G", "GGT": "G",
    "TCA": "S", "TCC": "S", "TCG": "S", "TCT": "S",
    "TTC": "F", "TTT": "F", "TTA": "L", "TTG": "L",
    "TAC": "Y", "TAT": "Y", "TAA": "_", "TAG": "_",
    "TGC": "C", "TGT": "C", "TGA": "_", "TGG": "W",
}  # fmt: skip


def translate_codon_naive(sequence: str) -> str:
    """Translate a nucleotide sequence, frame 0, no ORF/stop detection.

    Every 3 bases (including any past a stop codon) are translated; a
    trailing partial codon (sequence length not a multiple of 3) and any
    codon outside :data:`CODON_TABLE` (ambiguous bases, gaps) become ``X``.

    Args:
        sequence: nucleotide sequence (DNA, upper/lower-case accepted).

    Returns:
        Amino-acid string, same convention as bioconvert's ``fasta2faa``
        (stop codons as ``_``).

    Example::

        from sequana.codecs.fasta2faa import translate_codon_naive
        translate_codon_naive("ATGGCCTAA")
        'MA_'
    """
    seq = sequence.upper()
    return "".join(CODON_TABLE.get(seq[i : i + 3], "X") for i in range(0, len(seq), 3))


def fasta2faa(input_fasta: str, output_faa: str, width: int = 60) -> int:
    """Convert a nucleotide FASTA file to a protein FASTA (FAA) file.

    Args:
        input_fasta: path to a nucleotide FASTA file (``.fasta``/``.fa``,
            optionally ``.gz``  -- anything :class:`sequana.fasta.FastA` can
            read).
        output_faa: path to the output protein FASTA file.
        width: line-wrap width for the output sequence (default 60,
            matching bioconvert's convention).

    Returns:
        Number of sequences written.

    Example::

        from sequana.codecs.fasta2faa import fasta2faa
        fasta2faa("genes.fasta", "genes.faa")
    """
    written = 0
    with open(output_faa, "w") as fout:
        for record in FastA(input_fasta):
            name = record.name
            comment = getattr(record, "comment", "") or ""
            header = f"{name}\t{comment}" if comment else f"{name}\t"
            protein = translate_codon_naive(record.sequence)

            fout.write(f">{header}\n")
            fout.write("\n".join(textwrap.wrap(protein, width)) + "\n")
            written += 1

    return written
