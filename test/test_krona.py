from easydev import TempFile

from sequana import KronaMerger

from . import test_dir


def test_krona_merger():

    k1 = KronaMerger(f"{test_dir}/data/tsv/test_krona_k1.tsv")
    k2 = KronaMerger(f"{test_dir}/data/tsv/test_krona_k2.tsv")
    k1 += k2

    with TempFile(suffix=".tsv") as fh:
        df = k1.to_tsv(fh.name)
    assert all(df["count"] == [14043, 591, 184, 132])
    assert k1["Bacteria\tProteobacteria\tspecies1\n"] == 14043
