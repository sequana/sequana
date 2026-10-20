"""Parity of sequana_cov_rs.merge_rois against FilteredGenomeCov._merge_region."""
import numpy as np
import pandas as pd
import pytest

from sequana.bedtools import DoubleThresholds, FilteredGenomeCov

rs = pytest.importorskip("sequana_cov_rs")


def make_df(n=200_000, seed=0):
    rng = np.random.default_rng(seed)
    pos = np.arange(1, n + 1)
    z = rng.normal(0, 2, n)
    df = pd.DataFrame({"pos": pos, "cov": rng.poisson(30, n).astype(float), "rm": 30.0, "zscore": z})
    df = df[(df.zscore > 2) | (df.zscore < -2)].copy()  # same filter as get_rois
    df.index = df["pos"]
    df.insert(0, "chr", "chrTest")
    return df


@pytest.mark.parametrize("seed", [0, 1, 2])
def test_merge_rois_matches_python(seed):
    df = make_df(seed=seed)
    th = DoubleThresholds(-4, 4, 0.5, 0.5)
    py = pd.DataFrame(FilteredGenomeCov(df, th)._merge_region())
    out = rs.merge_rois(
        df["pos"].to_numpy(np.int64),
        df["cov"].to_numpy(),
        df["rm"].to_numpy(),
        df["zscore"].to_numpy(),
        th.low,
        th.high,
        1,
    )
    assert len(py) == len(out["start"]) > 0
    for col in ["start", "end", "size"]:
        np.testing.assert_array_equal(py[col].to_numpy(), out[col])
    for col in ["mean_cov", "mean_rm", "mean_zscore", "max_zscore", "max_cov", "log2_ratio"]:
        np.testing.assert_allclose(py[col].astype(float).to_numpy(), out[col], rtol=1e-12, equal_nan=True)


def _chrom(tmp_path_factory=None):
    from sequana import SequanaCoverage, sequana_data
    from sequana.bedtools_rs import ChromosomeCovRS

    gc = SequanaCoverage(sequana_data("virus.bed"))
    ch_py = gc[0]
    ch_rs = ChromosomeCovRS(gc, gc.chrom_names[0], gc.thresholds, gc.chunksize)
    for ch in (ch_py, ch_rs):
        ch.run(101, k=2, circular=False, force_models=False)
    return ch_py, ch_rs


def test_centralness_matches_python():
    import random

    random.seed(1)
    ch_py, ch_rs = _chrom()
    for t in (3, 4):
        assert ch_rs.get_centralness(t) == pytest.approx(ch_py.get_centralness(t), abs=1e-4)


def test_to_csv_matches_pandas():
    ch_py, ch_rs = _chrom()
    kw = dict(columns=["pos", "cov"], index=False, float_format="%.3g")
    assert ch_rs.to_csv(**kw) == ch_py.to_csv(**kw)
    assert ch_rs.to_csv(start=10, stop=500, **kw) == ch_py.to_csv(start=10, stop=500, **kw)


def test_fmt_g_random():
    rng = np.random.default_rng(0)
    values = np.concatenate([rng.normal(0, 1, 5000) * 10.0 ** rng.integers(-8, 8, 5000), [0.0, 999.5, 1e5, 0.0001]])
    df = pd.DataFrame({"x": values, "n": np.arange(len(values))})
    expected = df.to_csv(index=False, float_format="%.3g")
    got = rs.to_csv(["x", "n"], [values, df["n"].to_numpy(np.int64)], 3)
    assert got == expected


@pytest.mark.parametrize("seed", [0, 1])
def test_em_matches_python(seed):
    from sequana.mixture import EM
    from sequana.mixture_rs import EMRS

    rng = np.random.default_rng(seed)
    data = np.concatenate([rng.normal(1, 0.15, 70000), rng.normal(0.45, 0.1, 20000), rng.normal(2, 0.2, 10000)])
    a, b = EM(data), EMRS(data)
    a.estimate(k=2)
    b.estimate(k=2)
    for key in ("mus", "sigmas", "pis"):
        np.testing.assert_allclose(a.results[key], b.results[key], rtol=1e-8)
    for key in ("log_likelihood", "AIC", "BIC"):
        assert a.results[key] == pytest.approx(b.results[key], rel=1e-8)
    assert a.results.nfev == b.results.nfev and a.results.success == b.results.success
    np.testing.assert_allclose(a.gamma, b.gamma, rtol=1e-7, atol=1e-12)


@pytest.mark.parametrize("precision", [1, 2, 3, 6, 10])
def test_fmt_g_precisions(precision):
    rng = np.random.default_rng(precision)
    values = np.concatenate(
        [
            rng.normal(0, 1, 20000) * 10.0 ** rng.integers(-12, 12, 20000),
            rng.integers(0, 100000, 5000).astype(float),
            rng.integers(0, 1000, 5000) / 8.0,
            [
                0.0,
                -0.0,
                9.995,
                99.95,
                999.5,
                0.00009995,
                1e-5,
                1e15,
                np.inf,
                -np.inf,
                np.nan,
                5e-324,
                1.7976931348623157e308,
            ],
        ]
    )
    expected = pd.DataFrame({"x": values}).to_csv(index=False, float_format=f"%.{precision}g")
    assert rs.to_csv(["x"], [values], precision) == expected


@pytest.mark.parametrize("ncols, floats", [(3, False), (3, True), (4, True)])
def test_bed_reader_matches_pandas(tmp_path, ncols, floats):
    from sequana.bedtools import ChromosomeCov, SequanaCoverage
    from sequana.bedtools_rs import ChromosomeCovRS, SequanaCoverageRS

    rng = np.random.default_rng(0)
    lines = []
    for name, n in (("chrA", 3000), ("chrB", 1500), ("chrC", 800)):
        for pos in range(1, n + 1):
            cov = rng.poisson(30)
            fields = [name, str(pos), f"{cov:.2f}" if floats else str(cov)]
            if ncols == 4:
                fields.append(f"{rng.poisson(2):.1f}")
            lines.append("\t".join(fields))
    bed = tmp_path / "test.bed"
    bed.write_text("\n".join(lines) + "\n")

    a, b = SequanaCoverage(str(bed), chunksize=2000, quiet_progress=True), SequanaCoverageRS(str(bed), chunksize=2000)
    assert (a.chrom_names, a.positions, a.total_length) == (b.chrom_names, b.positions, b.total_length)
    for name in a.chrom_names:
        ca, cb = ChromosomeCov(a, name, a.thresholds, a.chunksize), ChromosomeCovRS(b, name, b.thresholds, b.chunksize)
        assert ca._mode == cb._mode
        while True:  # chrA is read in several chunks
            ca.next(), cb.next()
            if ca._df is None or cb._df is None:
                assert ca._df is None and cb._df is None
                break
            pd.testing.assert_frame_equal(ca._df, cb._df)
            if ca._mode == "memory":
                break
