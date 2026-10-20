import numpy as np
import pytest

from sequana.mixture import EM, norm_pdf


@pytest.mark.parametrize("mu, sigma", [(0, 1), (1.003, 0.066), (-5, 40.5), (0, 0), (0, -1), (0, np.nan)])
def test_norm_pdf_matches_scipy(mu, sigma):
    ss = pytest.importorskip("scipy.stats")
    x = np.concatenate([np.random.default_rng(0).normal(0, 3, 1000), [np.nan, np.inf, -np.inf, 0.0]])
    with np.errstate(all="ignore"):
        expected = ss.norm.pdf(x, mu, sigma)
    np.testing.assert_array_equal(norm_pdf(x, mu, sigma), expected)


def test_norm_pdf_scalar():
    assert norm_pdf(0.0) == pytest.approx(0.3989422804014327)
    assert np.isnan(norm_pdf(0.3, 1, 0))


def test_em_two_gaussians():
    rng = np.random.default_rng(1)
    em = EM(np.concatenate([rng.normal(1, 0.1, 7000), rng.normal(3, 0.2, 3000)]))
    em.estimate(k=2)
    assert sorted(round(m) for m in em.results.mus) == [1, 3]
