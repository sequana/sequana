import pytest

from sequana.genbank import GenBank

from . import test_dir


def test_features():

    gbk = GenBank(f"{test_dir}/data/genbank/JB409847.gbk")
    assert set(gbk.features) == set(["CDS", "source"])
    assert gbk.features  # repeats on purpose to test property

    gbk = GenBank(f"{test_dir}/data/genbank/unicycler.gbk")
    assert set(gbk.features) == set(["CDS", "source"])

    gbk = GenBank(f"{test_dir}/data/genbank/subsample.gb")
    assert set(gbk.features) == set(["CDS", "gene", "rRNA", "source", "tRNA"])


def test_features_parser():

    assert len(GenBank(f"{test_dir}/data/genbank/subsample.gb").genbank_features_parser()["NC_009641"]) == 220
    assert len(GenBank(f"{test_dir}/data/genbank/unicycler.gbk").genbank_features_parser()["1"]) == 40


def test_to_fasta(tmpdir):
    outname = tmpdir.join("test.fasta")
    gbk = GenBank(f"{test_dir}/data/genbank/JB409847.gbk")
    assert gbk.extract_fasta(f"{test_dir}/data/fasta/JB409847.fasta", features=["CDS"])
