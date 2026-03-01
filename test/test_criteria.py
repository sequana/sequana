import math

import pytest

from sequana.criteria import AIC, AICc, BIC


def test_AIC_with_likelihood():
    # AIC = 2k - 2*ln(L)
    result = AIC(1.0, 2, logL=False)
    expected = 2 * 2 - 2 * math.log(1.0)
    assert result == pytest.approx(expected)


def test_AIC_with_log_likelihood():
    # AIC = 2k + 2*logL (since L is already log)
    result = AIC(-10.0, 3, logL=True)
    expected = 2 * 3 + 2 * (-10.0)
    assert result == pytest.approx(expected)


def test_AICc():
    result = AICc(1.0, 2, 100, logL=False)
    aic = AIC(1.0, 2, logL=False)
    correction = 2 * 2 * (2 + 1.0) / (100 - 2 - 1.0)
    assert result == pytest.approx(aic + correction)


def test_BIC():
    result = BIC(1.0, 2, 100, logL=False)
    expected = -2 * math.log(1.0) + 2 * (math.log(100) - math.log(2 * math.pi))
    assert result == pytest.approx(expected)


def test_BIC_returns_float():
    result = BIC(1.0, 2, 100)
    assert isinstance(result, float)
