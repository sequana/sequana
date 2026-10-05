import pytest
from click.testing import CliRunner

from sequana.scripts.main import fasta_compare as script


@pytest.fixture
def fasta_files(tmp_path):
    seq1 = "ACGTACGTTTGACCATGCA" * 5
    seq2 = "GGGTTTCCCAAA" * 4
    seq3 = "ATATATGCGC" * 3
    ref = tmp_path / "ref.fa"
    ref.write_text(f">1\n{seq1}\n>2\n{seq2}\n>7\n{seq3}\n")
    asm = tmp_path / "asm.fa"
    asm.write_text(f">NC_0001.1\n{seq1}\n>NC_0002.1\n{seq2.lower()}\n>NC_9999.1\nTTTTTTTTTTGG\n")
    return ref, asm


def test_help():
    runner = CliRunner()
    results = runner.invoke(script.fasta_compare, ["--help"])
    assert results.exit_code == 0


def test_identical_files(fasta_files):
    ref, _ = fasta_files
    runner = CliRunner()
    results = runner.invoke(script.fasta_compare, [str(ref), str(ref)])
    assert results.exit_code == 0
    assert "1 == 1" in results.output


def test_different_names(fasta_files):
    ref, asm = fasta_files
    runner = CliRunner()
    results = runner.invoke(script.fasta_compare, [str(ref), str(asm)])
    # the two files differ (orphans on both sides), hence the exit code
    assert results.exit_code == 1
    assert "1 == NC_0001.1" in results.output
    assert "2 == NC_0002.1" in results.output
    assert "7" in results.output
    assert "NC_9999.1" in results.output


def test_strict_case(fasta_files):
    ref, asm = fasta_files
    runner = CliRunner()
    results = runner.invoke(script.fasta_compare, [str(ref), str(asm), "--strict-case"])
    assert results.exit_code == 1
    # the lower-case sequence is not matched anymore
    assert "2 == NC_0002.1" not in results.output


def test_rc_aware(tmp_path):
    from sequana.tools import reverse_complement

    seq = "ACGTACGTTTGACCATGCA" * 5
    f1 = tmp_path / "f1.fa"
    f1.write_text(f">1\n{seq}\n")
    f2 = tmp_path / "f2.fa"
    f2.write_text(f">NC_0001.1\n{reverse_complement(seq)}\n")

    runner = CliRunner()
    results = runner.invoke(script.fasta_compare, [str(f1), str(f2)])
    assert results.exit_code == 1

    results = runner.invoke(script.fasta_compare, [str(f1), str(f2), "--rc-aware"])
    assert results.exit_code == 0
    assert "1 == NC_0001.1" in results.output
    assert "reverse complement" in results.output


def test_output(fasta_files, tmp_path):
    ref, asm = fasta_files
    output = tmp_path / "report.tsv"
    runner = CliRunner()
    results = runner.invoke(script.fasta_compare, [str(ref), str(asm), "-o", str(output)])
    assert results.exit_code == 1
    lines = output.read_text().splitlines()
    assert lines[0].split("\t") == ["status", "names1", "names2", "length", "reverse_complement"]
    assert len([x for x in lines if x.startswith("match")]) == 2
    assert len([x for x in lines if x.startswith("orphan1")]) == 1
    assert len([x for x in lines if x.startswith("orphan2")]) == 1
