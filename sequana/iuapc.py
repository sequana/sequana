from sequana import logger
logger.name = __name__


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

