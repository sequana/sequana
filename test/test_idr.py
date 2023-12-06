import pytest

from sequana.idr import IDR

from . import test_dir


def test_idr(tmpdir):

    fout = tmpdir.join("test.png")

    for mode in ["narrow", "broad"]:
        idr = IDR(f"{test_dir}/data/IDR/{mode}.csv")
        assert idr.mode == mode
        len(idr)
        idr.plot_ranks(fout, savefig=True)
        idr.plot_scores(fout, savefig=True)
        idr.plot_rank_vs_idr_score(fout, savefig=True)
        idr.plot_idr_vs_peaks(fout, savefig=True)

        assert idr.IDR2score(0) == 1000
        assert idr.IDR2score(0.05) == 540
        assert idr.IDR2score(1) == 0
        assert idr.score2IDR(0) == 1
        assert pytest.approx(idr.score2IDR(540), 0.01) == 0.05
