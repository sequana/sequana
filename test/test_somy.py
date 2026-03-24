import numpy as np
import pandas as pd
import pytest

from sequana.somy import ClusterModels, SomyScore, plot_sliding_window_boxplot


def test_cluster_models_basic():
    mus = [10.0, 20.0, 30.0]
    sigmas = [1.0, 2.0, 3.0]
    pis = [0.5, 0.3, 0.2]
    cm = ClusterModels(mus, sigmas, pis)
    assert len(cm) == 3
    assert cm[0] == (10.0, 1.0, 0.5)
    assert cm[1] == (20.0, 2.0, 0.3)
    assert cm[2] == (30.0, 3.0, 0.2)


def test_cluster_models_get_std():
    mus = [10.0, 20.0]
    sigmas = [2.0, 4.0]
    pis = [0.6, 0.4]
    cm = ClusterModels(mus, sigmas, pis)
    expected = (0.6 * 2.0 + 0.4 * 4.0) / (0.6 + 0.4)
    assert abs(cm.get_std() - expected) < 1e-10


def test_cluster_models_cluster_well_separated():
    # Two gaussians far apart: only the dominant one (index 0) should be kept
    mus = [10.0, 100.0]
    sigmas = [1.0, 1.0]
    pis = [0.8, 0.2]
    cm = ClusterModels(mus, sigmas, pis)
    result = cm.cluster()
    assert isinstance(result, float)


def test_somy_score_init_no_file():
    ss = SomyScore(filename=None)
    assert ss.filename is None
    assert ss.window_size == 1000
    assert ss.df is None
    assert ss.save_em is True


def test_somy_score_init_custom_window():
    ss = SomyScore(filename=None, window_size=500)
    assert ss.window_size == 500


def _make_somy_df():
    """Build a minimal SomyScore._df for testing filter methods."""
    np.random.seed(42)
    data = {
        "chr": ["chr1"] * 40 + ["chr2"] * 40,
        "start": list(range(40)) + list(range(40)),
        "depth": np.concatenate(
            [
                np.random.normal(50, 5, 40),
                np.random.normal(100, 10, 40),
            ]
        ),
        "dist": list(range(40)) + list(range(40)),
        "tag": ["default"] * 80,
    }
    return pd.DataFrame(data)


def test_somy_score_remove_outliers():
    ss = SomyScore(filename=None)
    ss._df = _make_somy_df()
    n_before = len(ss._df)
    ss.remove_outliers(percentile=0.05)
    assert len(ss._df) < n_before


def test_somy_score_remove_flanks():
    ss = SomyScore(filename=None)
    ss._df = _make_somy_df()
    # window_size=1000, remove_flanking_regions_kb=10 -> N = ceil(10000/1000) = 10
    # rows with dist <= 10 should be removed
    n_before = len(ss._df)
    ss.remove_flanks(remove_flanking_regions_kb=10)
    assert len(ss._df) <= n_before
    assert (ss._df["dist"] > 10).all()


def test_somy_score_remove_low_depth():
    ss = SomyScore(filename=None)
    ss._df = _make_somy_df()
    ss.remove_low_depth(threshold=60)
    assert (ss._df["depth"] > 60).all()


def test_somy_score_boxplot_median_method():
    ss = SomyScore(filename=None)
    ss._df = _make_somy_df()
    # Should not raise; uses simple median normalisation
    ss.boxplot(method="median", normalise=True)


def test_somy_score_boxplot_mean_method():
    ss = SomyScore(filename=None)
    ss._df = _make_somy_df()
    ss.boxplot(method="mean", normalise=True)


def test_somy_score_boxplot_no_normalise():
    ss = SomyScore(filename=None)
    ss._df = _make_somy_df()
    ss.boxplot(normalise=False)


def test_somy_score_boxplot_muhat():
    ss = SomyScore(filename=None)
    ss._df = _make_somy_df()
    ss.boxplot(muhat=50.0, normalise=True)
    assert "error" in ss.info


def test_somy_score_somies_dataframe():
    ss = SomyScore(filename=None)
    ss._df = _make_somy_df()
    ss.boxplot(method="median", normalise=True)
    assert hasattr(ss, "somies")
    assert "estimated_somies" in ss.somies.columns
    assert "measured_somies" in ss.somies.columns
    assert set(ss.somies["estimated_somies"]).issubset({1, 2, 3, 4, 5})


def test_plot_sliding_window_boxplot():
    rng = np.random.default_rng(0)
    data = pd.Series(rng.normal(0, 1, 200))
    windows = plot_sliding_window_boxplot(data, window_size=20, step=1000, overlap_percentage=50)
    assert isinstance(windows, list)
    assert len(windows) > 0
    for w in windows:
        assert len(w) == 20


def test_plot_sliding_window_boxplot_too_large():
    data = pd.Series([1.0, 2.0, 3.0])
    with pytest.raises(ValueError, match="Window size too large"):
        plot_sliding_window_boxplot(data, window_size=100, step=1000, overlap_percentage=50)
