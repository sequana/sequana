import pytest

from sequana.cpg_islands import CpG, compute_cpg_content


def test_compute_cpg_content():
    result = compute_cpg_content("ACGTCG")
    assert "CpG count" in result
    assert "Observed/Expected CpG" in result
    assert "GC content" in result
    assert result["CpG count"] == 2


def test_compute_cpg_content_no_gc():
    # Sequence with no G or C
    result = compute_cpg_content("AAAA")
    assert result["CpG count"] == 0
    assert result["Observed/Expected CpG"] == 0
    assert result["GC content"] == 0.0


def test_cpg_calls_compute():
    result = CpG("ACGTCG")
    assert result["CpG count"] == 2


def test_compute_cpg_case_insensitive():
    lower_result = compute_cpg_content("acgtcg")
    upper_result = compute_cpg_content("ACGTCG")
    assert lower_result["CpG count"] == upper_result["CpG count"]
