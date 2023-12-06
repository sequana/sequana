from sequana.kmer import build_kmer, get_kmer


def test_build_kmer():
    res = build_kmer(length=3, letters="GC")
    assert len(res) == 8
    assert "CCG" in res
    assert "CCC" in res
    assert "CGC" in res
    assert "CGG" in res
    assert "GCC" in res
    assert "GCG" in res
    assert "GGC" in res
    assert "GGG" in res


def test_get_kmer():
    res = list(get_kmer("ACGTAAAA", k=4))
    assert res == ["ACGT", "CGTA", "GTAA", "TAAA", "AAAA"]
