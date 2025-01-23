from collections import Counter


def compute_melting_temperature_wallace_rule(sequence):
    """Compute mekting temperature Tm of a sequence using Wallace rule

    This rule is a quick estimation for short oligonucleotides (20-25 base pairs), based on [Marmu and Doty 1962]:

        Tm(Celsius) = 2 ({A} + {T}) + 4 ({G} + {C})

    Where A and T bases contribute 2°C each and G and C bases contribute 4°C each.

    This formula assumes standard conditions and is less accurate for longer sequences or those with unusual salt concentrations.
    """
    counter = Counter(sequence)
    return 2 * (counter["A"] + counter["T"]) + 4 * (counter["G"] + counter["C"])


def compute_melting_temperature_salt_adjusted(sequence):
    """compute melting temperature with salt adjustement

    This rule is a quick estimation for sequences greater than 14bp in length (Chester and Marshak 1993)

        Tm(Celsius) = 64.9 + 0.41 x %GC - 650 / (sequence length)

    This formula accounts for the stability conferred by GC content but does not account for secondary structures or mismatches.

    """
    counter = Counter(sequence)
    return 69.3 + 41 * (counter["G"] + counter["C"]) / len(sequence) - 650.0 / len(sequence)


# The DNA molecular weight assumes no modification of the terminal groups of the sequence.
molecular_weights_dna_bases = {"A": 313.21, "T": 304.2, "G": 329.21, "C": 289.18}

molecular_weights_rna_bases = {"A": 329.21, "U": 306.2, "G": 345.21, "C": 305.18}

# If the sequence is a single-stranded, synthesised oligonucleotide the values must be adjusted (removed phosphate group) molecular weight = molecular_weights_bases -61.96

# If the sequence is a single-stranded, cut by a restriction enzyme the value must be adjusted for the extra 5′-monophosphate left by most restriction enzymes by using: Molecular Weight = calculated molecular weight - 61.96 + 79.0


molecular_weights_amino_acids = {
    "A": 71.0788,
    "R": 156.1875,
    "N": 114.1038,
    "D": 115.0886,
    "C": 103.1388,
    "E": 129.1155,
    "Q": 128.1307,
    "G": 57.0519,
    "H": 137.1411,
    "I": 113.1594,
    "L": 113.1594,
    "K": 128.1741,
    "M": 131.1926,
    "F": 147.1766,
    "P": 97.1167,
    "S": 87.0782,
    "T": 101.1051,
    "W": 186.2132,
    "Y": 163.1760,
    "V": 99.1326,
    "U": 150.0388,
    "O": 237.3018,
}
