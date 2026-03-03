import pytest

from sequana.blast import blast_to_gff

from . import test_dir


def test_blast_to_gff(tmp_path):
    blast_input = tmp_path / "test.blastn"
    gff_output = tmp_path / "test.gff"

    # Write a minimal blast outfmt=6 file
    blast_input.write_text(
        "query1\tchrom1\t99.0\t100\t0\t0\t1\t100\t200\t300\t1e-50\t180\n"
        "query2\tchrom1\t98.0\t50\t1\t0\t1\t50\t600\t550\t1e-30\t95\n"  # reverse strand
    )

    blast_to_gff(str(blast_input), str(gff_output))

    content = gff_output.read_text()
    assert content.startswith("##gff-version 3\n")
    assert "chrom1" in content
    assert "query1" in content
    # forward strand
    assert "\t+\t" in content
    # reverse strand (sstart > send)
    assert "\t-\t" in content


def test_blast_to_gff_skip_comments(tmp_path):
    blast_input = tmp_path / "test.blastn"
    gff_output = tmp_path / "test.gff"

    blast_input.write_text(
        "# this is a comment\n"
        "\n"
        "query1\tchrom1\t99.0\t100\t0\t0\t1\t100\t200\t300\t1e-50\t180\n"
    )

    blast_to_gff(str(blast_input), str(gff_output))
    content = gff_output.read_text()
    lines = [l for l in content.splitlines() if not l.startswith("#")]
    assert len(lines) == 1
