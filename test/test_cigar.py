from sequana import Cigar, cigar


def test_cigar():

    cigar = "10S80M5I5I1D"  # the last 1D must be ignored, 5I5I should be understood
    c = Cigar(cigar)
    assert len(c) == 100

    assert c.as_tuple() == (("S", 10), ("M", 80), ("I", 5), ("I", 5), ("D", 1))
    assert c.as_sequence() == "S" * 10 + "M" * 80 + "I" * 10 + "D"
    assert c.as_dict() == {"S": 10, "M": 80, "I": 10, "D": 1}

    print(c)
    c.__repr__()

    c = Cigar("1S1S10M")
    c.compress()
    assert c.cigarstring == "2S10M"

    c = Cigar("1S10M")
    c.compress()
    assert c.cigarstring == "1S10M"

    c = Cigar("1S1S1S1S")
    c.compress()
    assert c.cigarstring == "4S"

    c = Cigar("1S2M1S")
    assert c.stats() == [2, 0, 0, 0, 2, 0, 0, 0, 0, 0]


def test_fetchers():
    cigar_list = [(0, 30), (1, 1), (2, 1), (3, 10), (0, 30), (4, 1000)]
    assert cigar.fetch_exon("chr1", 100, cigar_list) == [("chr1", 100, 130), ("chr1", 141, 171)]
    assert cigar.fetch_intron("chr1", 100, cigar_list) == [("chr1", 131, 141)]
    assert cigar.fetch_deletion("chr1", 100, cigar_list) == [("chr1", 130, 131)]
    assert cigar.fetch_clip("chr1", 100, cigar_list) == [("chr1", 171, 1171)]
    assert cigar.fetch_insertion("chr1", 100, cigar_list) == [("chr1", 130, 1)]

    # continue case
    assert cigar.fetch_exon("chr1", 100, [(40, 100)]) == []
    assert cigar.fetch_intron("chr1", 100, [(10, 10)]) == []

    """assert cigar.fetch_intron_string("chr1", 100, "30M1I1D30M30N30S") == \
        [['chr1', 162, 192]]
    assert cigar.fetch_deletion_string("chr1", 100, "30M1I1D30M30N30S") == \
        [['chr1', 131, 132]]
    assert cigar.fetch_insertion_string("chr1", 100, "30M1I1D30M30N30S") == \
        [['chr1', 130, 131]]
    assert cigar.fetch_exon_string("chr1", 100, "30M1I1D30M30N30S") == \
        [['chr1', 100, 130], ['chr1', 132, 162]]
    """
