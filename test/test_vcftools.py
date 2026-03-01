import pytest

from sequana.vcftools import compute_frequency, compute_strand_balance, strand_ratio

from . import test_dir


def test_strand_ratio():
    assert strand_ratio(3, 7) == pytest.approx(0.3)
    assert strand_ratio(5, 5) == pytest.approx(0.5)
    # Zero division should return 0
    assert strand_ratio(0, 0) == 0


def test_strand_ratio_above_half():
    # If division > 0.5, result should be 1 - division
    result = strand_ratio(8, 2)
    assert result == pytest.approx(0.2)


class _FakeRecord:
    """Minimal fake VCF record for testing."""

    def __init__(self, info):
        self.info = info


def test_compute_frequency_freebayes():
    record = _FakeRecord({"AO": [85], "DP": 100})
    result = compute_frequency(record)
    assert result == [pytest.approx(0.85)]


def test_compute_frequency_empty():
    # When neither AO nor VAF is present
    record = _FakeRecord({})
    result = compute_frequency(record)
    assert result == []


def test_compute_strand_balance():
    record = _FakeRecord({"SAF": [4], "SAR": [6]})
    result = compute_strand_balance(record)
    assert result == [pytest.approx(0.4)]


def test_compute_strand_balance_no_key():
    record = _FakeRecord({})
    result = compute_strand_balance(record)
    assert result == [0.5]
