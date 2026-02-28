from sequana.iuapc import *


def test_all():
    assert dna_bases == ("A", "C", "T", "G")
    assert dna_ambiguities["B"] == "[CGT]"
    assert dna_ambiguities_r["AG"] == "R"
    assert dna_complement["A"] == "T"
    assert codons["AUC"] == "I"
