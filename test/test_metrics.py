import pytest

from sequana.metrics import brukner_flexibility, compute_bendability, compute_helix_twist, helix_twist, lempel_ziv_complexity


def test_lempel_ziv_complexity_repetitive():
    # Highly repetitive sequence should have lower complexity
    low = lempel_ziv_complexity("AAAAAAAAAA")
    high = lempel_ziv_complexity("ACGTACGTGT")
    assert low < high


def test_lempel_ziv_complexity_returns_float():
    result = lempel_ziv_complexity("ACGT")
    assert isinstance(result, float)
    assert result > 0


def test_compute_bendability():
    seq = "ACGTA"
    scores = compute_bendability(seq, brukner_flexibility, window=3)
    assert len(scores) == len(seq) - 3 + 1
    # All trinucleotides in seq should be in the scale
    assert all(s is not None for s in scores)


def test_compute_bendability_unknown_trinuc():
    # Use a scale that doesn't cover all trinucleotides
    scores = compute_bendability("NNN", {}, window=3)
    assert scores == [None]


def test_compute_helix_twist():
    seq = "ACGT"
    scores = compute_helix_twist(seq, helix_twist)
    assert len(scores) == len(seq) - 1
    assert all(s is not None for s in scores)


def test_compute_helix_twist_unknown_dinuc():
    scores = compute_helix_twist("NN", {})
    assert scores == [None]
