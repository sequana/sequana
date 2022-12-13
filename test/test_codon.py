from sequana.codon import Codon


def test_find_start_position():

    c = Codon()

    # searches ATG on strand +
    assert c.find_start_codon_position("CCCATGAAA", 3, "+") == (3, "ATG")
    assert c.find_start_codon_position("CCCATGAAA", 0, "+") == (3, "ATG")

    # searches CAT on strand -
    # !! on strand -, returned position is shifted by 3 to start on the 3prime
    assert c.find_start_codon_position("AAACATGGG", 6, "-") == (3, "CAT")
    assert c.find_start_codon_position("AAACATGGG", 0, "-") == (3, "CAT")

    # searches stop TAG on strand +
    assert c.find_stop_codon_position("CCCTAG", 3, "+") == (3, "TAG")
    assert c.find_stop_codon_position("CCCTAG", 4, "+") == (3, "TAG")

    # searches stop TTA on strand -
    # !! on strand -, returned position is shifted by 3 to start on the 3prime
    assert c.find_stop_codon_position("CCCATTAAA", 6, "-") == (4, "TTA")
    assert c.find_stop_codon_position("CCCATTAAA", 0, "-") == (4, "TTA")
