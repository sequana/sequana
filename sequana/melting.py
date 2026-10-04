#
#  This file is part of Sequana software
#
#  Copyright (c) 2016-2026 - Sequana Development Team
#
#  Distributed under the terms of the 3-clause BSD license.
#  The full license is in the LICENSE file, distributed with this software.
#
#  website: https://github.com/sequana/sequana
#  documentation: http://sequana.readthedocs.io
#
##############################################################################
"""Melting temperature of nucleic acid duplexes with the nearest-neighbour model.

Nearest-neighbour parameters are provided for DNA/DNA duplexes (four sets) and for RNA/DNA hybrids. The RNA/DNA set
is what is needed to predict the stability of a DNA probe bound to an RNA target (see
:mod:`sequana.ribodesigner2`). The calculation follows the standard two-state nearest-neighbour model with the
entropy salt correction of SantaLucia (1998). Perfectly matched duplexes are supported (no mismatch, no dangling end).
The results agree with ``Bio.SeqUtils.MeltingTemp.Tm_NN`` (Biopython) for the same parameters; this is checked in the
tests.

========== ===================================================================================
Table      Reference
========== ===================================================================================
DNA_NN1    Breslauer, K.J. et al. PNAS 83, 3746-3750 (1986)
DNA_NN2    Sugimoto, N. et al. Nucleic Acids Res. 24, 4501-4505 (1996)
DNA_NN3    Allawi, H.T. and SantaLucia, J. Biochemistry 36, 10581-10594 (1997)
DNA_NN4    SantaLucia, J. and Hicks, D. Annu. Rev. Biophys. Biomol. Struct. 33, 415-440 (2004)
R_DNA_NN1  Sugimoto, N. et al. Biochemistry 34, 11211-11216 (1995), RNA/DNA hybrids
========== ===================================================================================

The salt correction is from SantaLucia, J. PNAS 95, 1460-1465 (1998).
"""
import math

from sequana.tools import reverse_complement

__all__ = ["TABLES", "tm_nn", "tm_probe_rna_dna", "tm_dna_dna"]

# Nearest-neighbour parameters: (delta H in kcal/mol, delta S in cal/(mol K)). A key "XY/ZW" is the dinucleotide XY
# of the reference strand (5' to 3') paired with its complement ZW (3' to 5'); the same step read from the other
# strand is found by reversing the key. The "init*" and "sym" keys are initiation and symmetry corrections.

DNA_NN1 = {  # Breslauer et al. (1986)
    "init": (0, 0),
    "init_A/T": (0, 0),
    "init_G/C": (0, 0),
    "init_oneG/C": (0, -16.8),
    "init_allA/T": (0, -20.1),
    "init_5T/A": (0, 0),
    "sym": (0, -1.3),
    "AA/TT": (-9.1, -24.0),
    "AT/TA": (-8.6, -23.9),
    "TA/AT": (-6.0, -16.9),
    "CA/GT": (-5.8, -12.9),
    "GT/CA": (-6.5, -17.3),
    "CT/GA": (-7.8, -20.8),
    "GA/CT": (-5.6, -13.5),
    "CG/GC": (-11.9, -27.8),
    "GC/CG": (-11.1, -26.7),
    "GG/CC": (-11.0, -26.6),
}

DNA_NN2 = {  # Sugimoto et al. (1996)
    "init": (0.6, -9.0),
    "init_A/T": (0, 0),
    "init_G/C": (0, 0),
    "init_oneG/C": (0, 0),
    "init_allA/T": (0, 0),
    "init_5T/A": (0, 0),
    "sym": (0, -1.4),
    "AA/TT": (-8.0, -21.9),
    "AT/TA": (-5.6, -15.2),
    "TA/AT": (-6.6, -18.4),
    "CA/GT": (-8.2, -21.0),
    "GT/CA": (-9.4, -25.5),
    "CT/GA": (-6.6, -16.4),
    "GA/CT": (-8.8, -23.5),
    "CG/GC": (-11.8, -29.0),
    "GC/CG": (-10.5, -26.4),
    "GG/CC": (-10.9, -28.4),
}

DNA_NN3 = {  # Allawi and SantaLucia (1997)
    "init": (0, 0),
    "init_A/T": (2.3, 4.1),
    "init_G/C": (0.1, -2.8),
    "init_oneG/C": (0, 0),
    "init_allA/T": (0, 0),
    "init_5T/A": (0, 0),
    "sym": (0, -1.4),
    "AA/TT": (-7.9, -22.2),
    "AT/TA": (-7.2, -20.4),
    "TA/AT": (-7.2, -21.3),
    "CA/GT": (-8.5, -22.7),
    "GT/CA": (-8.4, -22.4),
    "CT/GA": (-7.8, -21.0),
    "GA/CT": (-8.2, -22.2),
    "CG/GC": (-10.6, -27.2),
    "GC/CG": (-9.8, -24.4),
    "GG/CC": (-8.0, -19.9),
}

DNA_NN4 = {  # SantaLucia and Hicks (2004)
    "init": (0.2, -5.7),
    "init_A/T": (2.2, 6.9),
    "init_G/C": (0, 0),
    "init_oneG/C": (0, 0),
    "init_allA/T": (0, 0),
    "init_5T/A": (0, 0),
    "sym": (0, -1.4),
    "AA/TT": (-7.6, -21.3),
    "AT/TA": (-7.2, -20.4),
    "TA/AT": (-7.2, -21.3),
    "CA/GT": (-8.5, -22.7),
    "GT/CA": (-8.4, -22.4),
    "CT/GA": (-7.8, -21.0),
    "GA/CT": (-8.2, -22.2),
    "CG/GC": (-10.6, -27.2),
    "GC/CG": (-9.8, -24.4),
    "GG/CC": (-8.0, -19.9),
}

R_DNA_NN1 = {  # Sugimoto et al. (1995), the reference strand is the RNA
    "init": (1.9, -3.9),
    "init_A/T": (0, 0),
    "init_G/C": (0, 0),
    "init_oneG/C": (0, 0),
    "init_allA/T": (0, 0),
    "init_5T/A": (0, 0),
    "sym": (0, 0),
    "TT/AA": (-11.5, -36.4),
    "GT/CA": (-7.8, -21.6),
    "CT/GA": (-7.0, -19.7),
    "AT/TA": (-8.3, -23.9),
    "TG/AC": (-10.4, -28.4),
    "GG/CC": (-12.8, -31.9),
    "CG/GC": (-16.3, -47.1),
    "AG/TC": (-9.1, -23.5),
    "TC/AG": (-8.6, -22.9),
    "GC/CG": (-8.0, -17.1),
    "CC/GG": (-9.3, -23.2),
    "AC/TG": (-5.9, -12.3),
    "TA/AT": (-7.8, -23.2),
    "GA/CT": (-5.5, -13.5),
    "CA/GT": (-9.0, -26.1),
    "AA/TT": (-7.8, -21.9),
}

#: Parameter sets by name.
TABLES = {
    "DNA_NN1": DNA_NN1,
    "DNA_NN2": DNA_NN2,
    "DNA_NN3": DNA_NN3,
    "DNA_NN4": DNA_NN4,
    "R_DNA_NN1": R_DNA_NN1,
}

_COMPLEMENT = {"A": "T", "C": "G", "G": "C", "T": "A"}
_GAS_CONSTANT = 1.987  # cal / (K mol)


def _check(sequence):
    sequence = sequence.upper()
    if len(sequence) < 2:
        raise ValueError("The sequence must contain at least 2 bases.")
    if set(sequence) - set("ACGT"):
        raise ValueError("Only A, C, G and T are supported (no ambiguous bases).")
    return sequence


def tm_nn(target, na_mm=50.0, probe_nm=250.0, target_nm=0.0, table="R_DNA_NN1", self_complementary=None):
    """Melting temperature (Celsius) of a perfectly matched duplex, nearest-neighbour model.

    The duplex is made of ``target`` and its complement (the probe). For RNA/DNA hybrids the target is the RNA.

    :param target: sequence, 5' to 3', made of A, C, G and T (write U as T).
    :param na_mm: monovalent cation concentration in mM.
    :param probe_nm: concentration of the probe (higher concentrated strand) in nM.
    :param target_nm: concentration of the target (lower concentrated strand) in nM. 0 means a probe in large
        excess; the two strands are at equal concentration when ``target_nm == probe_nm``.
    :param table: name of a parameter set of :data:`TABLES` or a dictionary with the same structure.
    :param self_complementary: whether the sequence is self-complementary (palindromic). Detected from the
        sequence when None.
    :return: the melting temperature in Celsius.
    :raises ValueError: if the sequence is shorter than 2 bases or contains other characters than A, C, G, T.

    >>> from sequana.melting import tm_nn
    >>> round(tm_nn("GCGCGCGCGCGCGCGCGCGCGCGC", na_mm=100, probe_nm=1000), 1)
    80.0
    >>> round(tm_nn("CGTTGACGTAGCTAGCATCG", table="DNA_NN4", na_mm=50, probe_nm=250, target_nm=250), 1)
    56.3
    """
    if isinstance(table, str):
        try:
            table = TABLES[table]
        except KeyError:
            raise ValueError(f"Unknown table '{table}'. Choose from {sorted(TABLES)}.")
    target = _check(target)
    if self_complementary is None:
        self_complementary = target == reverse_complement(target)

    def parameter(key):
        return table.get(key, (0, 0))

    delta_h, delta_s = parameter("init")

    # duplex with no G/C pair, or at least one
    if "G" not in target and "C" not in target:
        extra = parameter("init_allA/T")
    else:
        extra = parameter("init_oneG/C")
    delta_h += extra[0]
    delta_s += extra[1]

    # penalty when the 5' end is T or the 3' end is A
    for penalty in (target.startswith("T"), target.endswith("A")):
        if penalty:
            delta_h += parameter("init_5T/A")[0]
            delta_s += parameter("init_5T/A")[1]

    # terminal base pairs
    ends = target[0] + target[-1]
    n_at = ends.count("A") + ends.count("T")
    n_gc = ends.count("G") + ends.count("C")
    delta_h += parameter("init_A/T")[0] * n_at + parameter("init_G/C")[0] * n_gc
    delta_s += parameter("init_A/T")[1] * n_at + parameter("init_G/C")[1] * n_gc

    # nearest neighbours
    for i in range(len(target) - 1):
        step = target[i : i + 2]
        key = step + "/" + "".join(_COMPLEMENT[b] for b in step)
        h, s = table[key] if key in table else table[key[::-1]]
        delta_h += h
        delta_s += s

    if self_complementary:
        k = probe_nm * 1e-9
        delta_h += parameter("sym")[0]
        delta_s += parameter("sym")[1]
    else:
        k = (probe_nm - target_nm / 2.0) * 1e-9

    # salt correction on the entropy (SantaLucia 1998)
    delta_s += 0.368 * (len(target) - 1) * math.log(na_mm * 1e-3)

    return 1000.0 * delta_h / (delta_s + _GAS_CONSTANT * math.log(k)) - 273.15


def tm_probe_rna_dna(probe, na_mm=100.0, probe_nm=1000.0):
    """Melting temperature (Celsius) of a DNA probe bound to its complementary RNA target.

    Uses the RNA/DNA hybrid parameters of Sugimoto et al. (1995). The probe is assumed to be in large excess
    over the target.

    :param probe: DNA probe sequence, 5' to 3'. The RNA target is its reverse complement.
    :param na_mm: monovalent cation concentration in mM.
    :param probe_nm: concentration of the probe in nM.
    :raises ValueError: if the probe contains other characters than A, C, G, T.

    >>> from sequana.melting import tm_probe_rna_dna
    >>> round(tm_probe_rna_dna("AGCCCTCATTCTTTACGAAGCCCAGCGTCAGCGG"), 1)
    76.4
    """
    probe = _check(probe)
    return tm_nn(reverse_complement(probe), na_mm=na_mm, probe_nm=probe_nm, table="R_DNA_NN1")


def tm_dna_dna(sequence, na_mm=50.0, strand_nm=250.0, table="DNA_NN4"):
    """Melting temperature (Celsius) of a DNA oligonucleotide annealed to its perfect complement.

    Both strands are assumed to be at the same concentration (for instance a primer and its template).

    :param sequence: DNA sequence, 5' to 3'.
    :param na_mm: monovalent cation concentration in mM.
    :param strand_nm: concentration of each strand in nM.
    :param table: DNA/DNA parameter set, ``"DNA_NN1"`` to ``"DNA_NN4"`` (default: SantaLucia and Hicks, 2004).
    :raises ValueError: if the sequence contains other characters than A, C, G, T.

    >>> from sequana.melting import tm_dna_dna
    >>> round(tm_dna_dna("CGTTGACGTAGCTAGCATCG"), 1)
    56.3
    """
    if table == "R_DNA_NN1":
        raise ValueError("Use tm_probe_rna_dna for RNA/DNA hybrids.")
    return tm_nn(_check(sequence), na_mm=na_mm, probe_nm=strand_nm, target_nm=strand_nm, table=table)
