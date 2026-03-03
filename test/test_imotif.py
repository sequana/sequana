import pytest

from sequana.imotif import IMotif

from . import test_dir

fasta_file = f"{test_dir}/data/fasta/measles.fa"


def test_imotif_run():
    im = IMotif(fasta_file, min_tract=3, max_loop=7)
    im.run()
    assert "seqid" in im.df.columns
    assert "sequence" in im.df.columns


def test_imotif_to_bed(tmp_path):
    im = IMotif(fasta_file, min_tract=3, max_loop=7)
    im.run()
    if not im.df.empty:
        bed_out = tmp_path / "imotif.bed"
        im.to_bed(str(bed_out))
        assert bed_out.exists()


def test_imotif_to_bed_empty_raises(tmp_path):
    im = IMotif(fasta_file)
    with pytest.raises(ValueError):
        im.to_bed(str(tmp_path / "imotif_empty.bed"))
