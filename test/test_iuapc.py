from sequana.iuapc import *


def test_all():
    assert dna_bases == ("A", "C", "T", "G")
    assert dna_ambiguities["B"] == "[CGT]"
    assert dna_ambiguities_r["AG"] == "R"
    assert dna_complement["A"] == "T"
    assert codons["AUC"] == "I"
    assert amino_acids["A"] == ("Ala", "Alanine")
    assert amino_acids["M"] == ("Met", "Methionine")
    assert len(amino_acids) == 20
    assert exotic_amino_acids["X"] == ("Xaa", "Any residue")
    assert exotic_amino_acids["U"] == ("Sec", "Selenocysteine")
    assert exotic_amino_acids["O"] == ("Pyl", "Pyrrolysine")
    assert dna_bases_names["A"] == "Adenine"
    assert dna_bases_names["N"] == "Unknown"
    assert len(dna_bases) == 4
    assert len(codons) == 64
