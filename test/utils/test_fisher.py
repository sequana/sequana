import pytest

from sequana.utils.fisher import fisher_exact

# values were computed with scipy
# scipy.stats.fisher_exact(table=[[45, 40], [10, 5]], alternative="greater")
# scipy.stats.fisher_exact(table=[[45, 40], [10, 5]], alternative="less")
# scipy.stats.fisher_exact(table=[[45, 40], [10, 5]], alternative="two-sided")


def test_right():
    assert fisher_exact([[45, 40], [10, 5]], "greater") == pytest.approx(0.8985632205665318, 1e-9)


def test_left():
    assert fisher_exact([[45, 40], [10, 5]], "less") == pytest.approx(0.24249234232734680, 1e-9)


def test_both():
    assert fisher_exact([[45, 40], [10, 5]], "two-sided") == pytest.approx(0.4048121346770015, 1e-9)


def test_others():
    # strong asymetry
    fisher_exact([[45, 40], [10, 500000]], "two-sided")
    fisher_exact([[45, 40], [10, 500000]], "less")
    fisher_exact([[45, 40], [10, 500000]], "greater")
