import pytest

from sequana.wigtools import yield_wig_by_chromosome  # replace with actual import

# Fixtures for inline test files
wig_variable = """\
variableStep chrom=chr1
100 1.0
101 2.0
102 3.0
variableStep chrom=chr2
100 1.0
101 2.0
"""

wig_fixed = """\
fixedStep chrom=chr2 start=200 step=2
4.0
5.0
6.0
fixedStep chrom=chr3 start=200 step=2
4.0
5.0
6.0
"""


def test_yield_wig_by_chromosome(tmp_path):
    # Write variableStep test file
    file_var = tmp_path / "var.wig"
    file_var.write_text(wig_variable)

    result = list(yield_wig_by_chromosome(file_var))
    assert len(result) == 2
    chrom, entries = result[0]
    assert chrom == "chr1"
    assert entries == [(100, 1.0), (101, 2.0), (102, 3.0)]

    # Write fixedStep test file
    file_fix = tmp_path / "fix.wig"
    file_fix.write_text(wig_fixed)

    result = list(yield_wig_by_chromosome(file_fix))
    assert len(result) == 2
    chrom, entries = result[0]
    assert chrom == "chr2"
    assert entries == [(200, 4.0), (202, 5.0), (204, 6.0)]
