import pytest

from sequana.palindromes import Palindromes

from . import test_dir

fasta_file = f"{test_dir}/data/fasta/measles.fa"


def test_palindromes_run(tmp_path):
    p = Palindromes(fasta_file, min_len=4, max_len=6)
    p.run()
    # dataframe should have expected columns
    assert "seqid" in p.df.columns
    assert "sequence" in p.df.columns


def test_palindromes_is_palindrome():
    p = Palindromes(fasta_file)
    assert p.is_palindrome("ATAT")  # AT paired with AT reverse complement
    assert not p.is_palindrome("AAAA")


def test_palindromes_to_bed(tmp_path):
    p = Palindromes(fasta_file, min_len=4, max_len=6)
    p.run()
    if not p.df.empty:
        bed_out = tmp_path / "palindromes.bed"
        p.to_bed(str(bed_out))
        assert bed_out.exists()


def test_palindromes_to_bed_empty_raises(tmp_path):
    p = Palindromes(fasta_file)
    with pytest.raises(ValueError):
        p.to_bed(str(tmp_path / "palindromes_empty.bed"))
