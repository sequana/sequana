import pytest

from sequana.melting import TABLES, tm_dna_dna, tm_nn, tm_probe_rna_dna
from sequana.tools import reverse_complement

# Reference values computed with Bio.SeqUtils.MeltingTemp.Tm_NN(seq, nn_table=R_DNA_NN1, Na=100, dnac1=1000, dnac2=0)
# (target = RNA sequence).
REFERENCE = {
    "AGCCCTCATTCTTTACGAAGCCCAGCGTCAGCGGCAAAATGCAGGCATGTACCTGCG": 78.56352472473304,
    "GCGCGCGCGCGCGCGCGCGCGCGC": 80.0158660203279,
    "ATATATATATATATATATAT": 32.351596454552464,
    "TTTTTAAAAACCCCCGGGGG": 58.726358739463535,
}


@pytest.mark.parametrize("target,expected", REFERENCE.items())
def test_tm_nn_reference(target, expected):
    assert tm_nn(target, na_mm=100, probe_nm=1000) == pytest.approx(expected, abs=1e-6)


def test_tm_nn_matches_biopython_on_random_sequences():
    mt = pytest.importorskip("Bio.SeqUtils.MeltingTemp")
    import random

    rng = random.Random(1)
    for _ in range(200):
        seq = "".join(rng.choice("ACGT") for _ in range(rng.randint(20, 80)))
        for na, conc in [(50, 250), (100, 1000), (300, 500)]:
            expected = mt.Tm_NN(seq, nn_table=mt.R_DNA_NN1, Na=na, dnac1=conc, dnac2=0)
            assert tm_nn(seq, na_mm=na, probe_nm=conc) == pytest.approx(expected, abs=1e-6)


def test_tm_probe_rna_dna_uses_reverse_complement():
    probe = "AGCCCTCATTCTTTACGAAGCCCAGCGTCAGCGG"
    assert tm_probe_rna_dna(probe, 100, 1000) == pytest.approx(tm_nn(reverse_complement(probe), 100, 1000))


def test_tm_depends_on_salt_and_gc():
    assert tm_nn("GC" * 20, 100, 1000) > tm_nn("AT" * 20, 100, 1000)
    assert tm_nn("ACGT" * 10, 300, 1000) > tm_nn("ACGT" * 10, 50, 1000)


def test_invalid_sequences():
    with pytest.raises(ValueError):
        tm_nn("ACGTN")
    with pytest.raises(ValueError):
        tm_nn("A")
    with pytest.raises(ValueError):
        tm_probe_rna_dna("ACGUN")


# Bio.SeqUtils.MeltingTemp.Tm_NN(seq, nn_table=table, Na=50, dnac1=250, dnac2=250, selfcomp=<palindrome>)
DNA_SEQUENCES = [
    "CGTTGACGTAGCTAGCATCG",
    "GAATTCGAATTCGAATTC",
    "ATATATATATAT",
    "TTTTTTTTTTGCGCGCGC",
    "AAAAAAAAAAAAAAAAAAAT",
]
DNA_REFERENCE = {
    "DNA_NN1": [63.8253415181, 55.4559171359, 6.6540035605, 71.9580221649, 53.038080084],
    "DNA_NN2": [62.8167534991, 49.8105513673, 6.5907778991, 60.8750908392, 43.2456010105],
    "DNA_NN3": [56.382716006, 44.6964714367, 8.9633581612, 55.6539072193, 38.6975605516],
    "DNA_NN4": [56.3089248606, 44.4297097933, 8.8567122355, 55.4967694709, 37.8514441437],
}


@pytest.mark.parametrize("table", sorted(DNA_REFERENCE))
def test_tm_dna_dna_reference(table):
    # the 2nd and 3rd sequences are palindromes: the symmetry correction and single-strand concentration apply
    values = [tm_dna_dna(seq, na_mm=50, strand_nm=250, table=table) for seq in DNA_SEQUENCES]
    assert values == pytest.approx(DNA_REFERENCE[table], abs=1e-6)


@pytest.mark.parametrize("table", ["DNA_NN1", "DNA_NN2", "DNA_NN3", "DNA_NN4"])
def test_dna_dna_matches_biopython_on_random_sequences(table):
    mt = pytest.importorskip("Bio.SeqUtils.MeltingTemp")
    import random

    rng = random.Random(2)
    for _ in range(150):
        seq = "".join(rng.choice("ACGT") for _ in range(rng.randint(10, 40)))
        selfcomp = seq == reverse_complement(seq)
        expected = mt.Tm_NN(seq, nn_table=getattr(mt, table), Na=80, dnac1=500, dnac2=100, selfcomp=selfcomp)
        assert tm_nn(seq, na_mm=80, probe_nm=500, target_nm=100, table=table) == pytest.approx(expected, abs=1e-6)


def test_self_complementary_detection():
    palindrome = "GAATTC"
    assert tm_nn(palindrome, table="DNA_NN4", self_complementary=True) != tm_nn(
        palindrome, table="DNA_NN4", self_complementary=False
    )
    assert tm_nn(palindrome, table="DNA_NN4") == tm_nn(palindrome, table="DNA_NN4", self_complementary=True)


def test_unknown_table_and_rna_dna_guard():
    with pytest.raises(ValueError):
        tm_nn("ACGTACGT", table="nope")
    with pytest.raises(ValueError):
        tm_dna_dna("ACGTACGT", table="R_DNA_NN1")


def test_custom_table_dictionary():
    table = dict(TABLES["DNA_NN4"])
    assert tm_nn("ACGTACGTAC", table=table) == pytest.approx(tm_nn("ACGTACGTAC", table="DNA_NN4"))


def test_tables_complete():
    # 10 unique dinucleotide steps for DNA/DNA, 16 for RNA/DNA, plus 7 initiation/symmetry keys
    assert all(len(TABLES[name]) == 17 for name in ["DNA_NN1", "DNA_NN2", "DNA_NN3", "DNA_NN4"])
    assert len(TABLES["R_DNA_NN1"]) == 23
