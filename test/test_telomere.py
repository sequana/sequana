import os

from sequana.telomere import Telomere

from . import test_dir


def test_telomere(tmpdir):
    path = tmpdir.mkdir("temp")

    # 2 telomeres
    filename = f"{test_dir}/data/fasta/test_telomere.fasta"
    t = Telomere(filename)
    candidates = t.find_representative_kmers(t.fasta.sequences[0])
    assert "TAACCC" in [x[0] for x in candidates]
    t.run(tag=None)

    # no telomeres
    filename = f"{test_dir}/data/fasta/measles.fa"
    t = Telomere(filename)
    candidates = t.find_representative_kmers(t.fasta.sequences[0])
    assert len(candidates) == 0
    t.run(tag=None)
