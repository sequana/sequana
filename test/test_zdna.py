import pytest

from sequana.zdna import ZDNA

from . import test_dir

fasta_file = f"{test_dir}/data/fasta/measles.fa"


def test_zdna_run():
    z = ZDNA(fasta_file, motif="CG", min_repeats=2)
    z.run()
    assert "seqid" in z.df.columns
    assert "sequence" in z.df.columns


def test_zdna_to_bed(tmp_path):
    z = ZDNA(fasta_file, motif="CG", min_repeats=2)
    z.run()
    if not z.df.empty:
        bed_out = tmp_path / "zdna.bed"
        z.to_bed(str(bed_out))
        assert bed_out.exists()


def test_zdna_to_bed_empty_raises(tmp_path):
    z = ZDNA(fasta_file)
    with pytest.raises(ValueError):
        z.to_bed(str(tmp_path / "zdna_empty.bed"))
