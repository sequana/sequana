def lempel_ziv_complexity(seq):
    """

    Lempel-Ziv complexity is a measure of sequence compressibility—how "random" or "repetitive" a sequence is. In the context of DNA, it captures how diverse the patterns are. A low LZ complexity suggests repeats or low complexity regions, while a high LZ complexity suggests a more random or information-rich sequence.



    """

    seq = seq.upper()
    n = len(seq)
    i = 0
    complexity = 1
    dictionary = set()
    w = seq[0]

    while i < n - 1:
        i += 1
        w += seq[i]
        if w not in dictionary:
            dictionary.add(w)
            complexity += 1
            w = ""
            if i + 1 < n:
                w = seq[i + 1]
                i += 1
    return complexity / len(seq)


def compute_bendability(seq, scale, window=3):
    """

    DNA bendability measures the local flexibility of a DNA sequence—its ability to bend or deform. This is a crucial feature in:


            Brukner et al. (1995)

    """

    seq = seq.upper()
    scores = []
    for i in range(len(seq) - window + 1):
        tri = seq[i : i + window]
        score = scale.get(tri, None)
        if score is not None:
            scores.append(score)
        else:
            scores.append(None)
    return scores


brukner_flexibility = {
    "AAA": 0.069,
    "AAC": 0.055,
    "AAG": 0.059,
    "AAT": 0.063,
    "ACA": 0.072,
    "ACC": 0.054,
    "ACG": 0.051,
    "ACT": 0.058,
    "AGA": 0.070,
    "AGC": 0.056,
    "AGG": 0.058,
    "AGT": 0.064,
    "ATA": 0.066,
    "ATC": 0.053,
    "ATG": 0.060,
    "ATT": 0.065,
    "CAA": 0.073,
    "CAC": 0.057,
    "CAG": 0.061,
    "CAT": 0.067,
    "CCA": 0.070,
    "CCC": 0.052,
    "CCG": 0.050,
    "CCT": 0.057,
    "CGA": 0.069,
    "CGC": 0.055,
    "CGG": 0.056,
    "CGT": 0.062,
    "CTA": 0.065,
    "CTC": 0.051,
    "CTG": 0.059,
    "CTT": 0.066,
    "GAA": 0.068,
    "GAC": 0.054,
    "GAG": 0.058,
    "GAT": 0.064,
    "GCA": 0.071,
    "GCC": 0.053,
    "GCG": 0.050,
    "GCT": 0.057,
    "GGA": 0.067,
    "GGC": 0.055,
    "GGG": 0.057,
    "GGT": 0.063,
    "GTA": 0.064,
    "GTC": 0.052,
    "GTG": 0.059,
    "GTT": 0.065,
    "TAA": 0.067,
    "TAC": 0.053,
    "TAG": 0.060,
    "TAT": 0.066,
    "TCA": 0.069,
    "TCC": 0.054,
    "TCG": 0.052,
    "TCT": 0.058,
    "TGA": 0.068,
    "TGC": 0.056,
    "TGG": 0.059,
    "TGT": 0.064,
    "TTA": 0.065,
    "TTC": 0.051,
    "TTG": 0.060,
    "TTT": 0.067,
}

helix_twist = {
    "AA": 36.0,
    "AC": 34.5,
    "AG": 36.6,
    "AT": 32.2,
    "CA": 35.6,
    "CC": 33.1,
    "CG": 30.0,
    "CT": 35.2,
    "GA": 35.7,
    "GC": 36.9,
    "GG": 32.8,
    "GT": 34.3,
    "TA": 34.0,
    "TC": 35.4,
    "TG": 36.1,
    "TT": 36.0,
}


def compute_helix_twist(seq, scale):
    """

    Helix twist refers to the rotational angle between adjacent base pairs in the DNA double helix, typically measured in degrees. It reflects how tightly the DNA is twisted and is important for:

    - DNA supercoiling
    - Protein-DNA interactions
    - Nucleosome positioning
    - Local helical structure variations (e.g., bends, kinks)

    Typical Values: Canonical B-DNA: ~36° per base pair step
    Varies with sequence: ~32°–40° depending on dinucleotide
    """

    seq = seq.upper()
    scores = []
    for i in range(len(seq) - 1):
        dinuc = seq[i : i + 2]
        score = scale.get(dinuc, None)
        if score is not None:
            scores.append(score)
        else:
            scores.append(None)
    return scores
