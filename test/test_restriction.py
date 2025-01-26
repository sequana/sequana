from sequana.restriction import *


def test_restriction():
    # Find restriction sites
    dna_sequence = "GAATTCTAGCGGCCGCGAATTCGGTACC"
    restriction_sites = find_restriction_sites(dna_sequence, restriction_enzymes)

    assert restriction_sites["EcoRI"] == [1, 17]
    assert restriction_sites["BamHI"] == []
    assert restriction_sites["HindIII"] == []
    assert restriction_sites["NotI"] == [9]
