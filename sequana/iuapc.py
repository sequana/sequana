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
import colorlog
logger = colorlog.getLogger(__name__)



"""
There are 20 amino acids (Stored in :ref:`amino_acids`) such as alanine.
Each of them has a name, a 3-letter codes and a one-letter codes (the IUAPC
code, after the International Union of Pure and Applied chemistry committee).

`ref:`amino_acids is a dictionary with keys as the one-letter code. The values
are tuples with the 3-letter code and full name.

amino_acids['A'] returns ('Ala', 'Alanine')

Som a protein can be written as a succession of letters made of the keys of the
dictionary. It is then easy to check the validity of a protein sequence.



"""

# DNA bases
dna_bases = ("A", "C", "T", "G")

# DNA bases names

dna_bases_names = {
    "A": "Adenine",
    "T": "Thymidine",
    "U": "Uridine",
    "G": "Guanidine",
    "C": "Cytidine",
    "Y": "pYrimidine",
    "R": "puRine",
    "S": "Strong",
    "W": "Weak",
    "K": "Keto",
    "M": "aMino",
    "B": "not A",
    "D": "not C",
    "H": "not G",
    "V": "not T/U",
    "N": "Unknown"}

# DNA bases represented
dna_ambiguities = {
    "A": "A", 
    "C": "C",
    "G": "G",
    "T": "T",
    "R": "[GA]",
    "Y": "[CT]",
    "M": "[AC]",
    "K": "[GT]",
    "S": "[GC]",
    "W": "[AT]",
    "N": "[ACGT]",
    "B": "[CGT]",
    "D": "[AGT]",
    "H": "[ACT]",
    "V": "[ACG]"
}

# IUPAC degeneracies. Complementary bases
dna_complement = {
    'A': 'T',
    'B': 'V',
    'C': 'G',
    'D': 'H',
    'G': 'C',
    'H': 'D',
    'K': 'M',
    'M': 'K',
    'N': 'N',
    'R': 'Y',
    'S': 'S',
    'T': 'A',
    'V': 'B',
    'W': 'W',
    'X': 'X',
    'Y': 'R'}

codons = {
    "UUU":"F", "UUC":"F","UUA":"L", "UUG":"L",
    "CUU":"L", "CUC":"L","CUA":"L", "CUG":"L",
    "AUU":"I", "AUC":"I","AUA":"I", "AUG":"M",
    "GUU":"V", "GUC":"V","GUA":"V", "GUG":"V",
    "UCU":"S", "UCC":"S","UCA":"S", "UCG":"S",
    "CCU":"P", "ACC":"P","CCA":"P", "CCG":"P",
    "ACU":"T", "ACC":"T","ACA":"T", "ACG":"T",
    "GCU":"A", "GCC":"A","GCA":"A", "GCG":"A",
    "UAU":"Y", "UAC":"Y","UAA":"*", "UAG":"*",
    "CAU":"H", "CAC":"H","CAA":"Q", "CAG":"Q",
    "AAU":"N", "AAC":"N","AAA":"K", "AAG":"K",
    "GAU":"D", "GAC":"D","GAA":"E", "GAG":"E",
    "UGU":"C", "UGC":"C","UGA":"*", "UGG":"W",
    "CGU":"R", "CGC":"R","CGA":"R", "CGG":"R",
    "AGU":"S", "AGC":"S","AGA":"R", "AGG":"R",
    "GGU":"G", "GGC":"G","GGA":"G", "GGG":"G",
    }


amino_acids = {
"A": ('Ala', 'Alanine'),
"R": ('Arg', 'Arginine'),
"N": ('Asn', 'Asparagine'),
"D": ('Asp', 'Aspartic acid'),
"C": ('Cys', 'Cysteine'),
"Q": ('Gln', 'Glutamine'),
"E": ('Glu', 'Glutamic acid'),
"G": ('Gly', 'Glycine'),
"H": ('His', 'Histidine'),
"I": ('lle', 'Isoleucine'),
"L": ('Leu', 'Leucine'),
"K": ('Lys', 'Lysine'),
"M": ('Met', 'Methionine'),
"F": ('Phe', 'Pheline'),
"P": ('Pro', 'Proline'),
"S": ('Ser', 'Serine'),
"T": ('Thr', 'Threonine'),
"W": ('Trp', 'Tryptophan'),
"Y": ('Tyr', 'Tyrosine'),
"V": ('Val', 'Valine')
}

# B and Z codes indicated ambiguous amino acd
# Pyl and Sec are specified by the UAG  (Pyl) and UGA (Sec) stop codons in a
# specific context
exotic_amino_acids = {
"B": ("Asn or Asp", "Asparagine or aspartic acis"),
"J": ("Xle", "Isoleucine or leucine"),
"O": ("Pyl", "Pyrrolysine"),
"U": ("Sec", "Selenocysteine"),
"Z": ("Gln or Glu", "Glutamine or glutamic acid"),
"X": ("Xaa", "Any residue"),
"--": ("gap", "gap"),


}









