# extract from https://github.com/painyeph/FishersExactTest/tree/master
# we are only interested in the 2-tail test that is 150 times faster than scipy.

# changes made.
# - only a subset of essential functions were copied (test1l, test1r, test1t)
# - The 'all' flavor was dropped
# - created a single entry point function similar to scipy syntax.
# - uses python caching instead of own caching (same performances but more reliable on long term)
from math import exp, lgamma, log

try:
    from functools import cache
except ImportError:
    from functools import lru_cache as cache


@cache
def lgamma(z):
    x = 0.1659470187408462e-06 / (z + 7)
    x += 0.9934937113930748e-05 / (z + 6)
    x -= 0.1385710331296526 / (z + 5)
    x += 12.50734324009056 / (z + 4)
    x -= 176.6150291498386 / (z + 3)
    x += 771.3234287757674 / (z + 2)
    x -= 1259.139216722289 / (z + 1)
    x += 676.5203681218835 / z
    x += 0.9999999999995183
    return log(x) - 5.58106146679532777 - z + (z - 0.5) * log(z + 6.5)


def _maxn():
    l = 1
    n = 2
    h = float("inf")
    while l < n:
        if abs(lgamma(n + 1) - lgamma(n) - log(n)) >= 1:
            h = n
        else:
            l = n
        n = (l + min(h, l * 3)) // 2
    return n


MAXN = _maxn()


def fisher_exact(array_2by2, alternative):
    """Fast Fisher test alternative to scipy"""
    a = array_2by2[0][0]
    b = array_2by2[0][1]
    c = array_2by2[1][0]
    d = array_2by2[1][1]
    if alternative == "less":
        return test_left(a, b, c, d)
    elif alternative == "greater":
        return test_right(a, b, c, d)
    elif alternative == "two-sided":
        return test_both(a, b, c, d)
    else:  # pragma: no cover
        raise ValueError("alternative must be less, greater or two-sided")


# two-tails
def test_both(a, b, c, d):
    return exp(-mlnTest2t(a, a + b, a + c, a + b + c + d))


def mlnTest2t(a, ab, ac, abcd):
    if 0 > a or a > ab or a > ac or ab + ac > abcd + a:
        raise ValueError("invalid contingency table")
    if abcd > MAXN:
        raise OverflowError("the grand total of contingency table is too large")
    a_min = max(0, ab + ac - abcd)
    a_max = min(ab, ac)
    if a_min == a_max:
        return 0.0
    p0 = lgamma(ab + 1) + lgamma(ac + 1) + lgamma(abcd - ac + 1) + lgamma(abcd - ab + 1) - lgamma(abcd + 1)
    pa = lgamma(a + 1) + lgamma(ab - a + 1) + lgamma(ac - a + 1) + lgamma(abcd - ab - ac + a + 1)
    st = 1.0
    if ab * ac < a * abcd:
        for i in range(min(a - 1, int(round(ab * ac / abcd))), a_min - 1, -1):
            pi = lgamma(i + 1) + lgamma(ab - i + 1) + lgamma(ac - i + 1) + lgamma(abcd - ab - ac + i + 1)
            if pi < pa:
                continue
            st_new = st + exp(pa - pi)
            if st_new == st:
                break
            st = st_new
        for i in range(a + 1, a_max + 1):
            pi = lgamma(i + 1) + lgamma(ab - i + 1) + lgamma(ac - i + 1) + lgamma(abcd - ab - ac + i + 1)
            st_new = st + exp(pa - pi)
            if st_new == st:
                break
            st = st_new
    else:
        for i in range(a - 1, a_min - 1, -1):
            pi = lgamma(i + 1) + lgamma(ab - i + 1) + lgamma(ac - i + 1) + lgamma(abcd - ab - ac + i + 1)
            st_new = st + exp(pa - pi)
            if st_new == st:
                break
            st = st_new
        for i in range(max(a + 1, int(round(ab * ac / abcd))), a_max + 1):
            pi = lgamma(i + 1) + lgamma(ab - i + 1) + lgamma(ac - i + 1) + lgamma(abcd - ab - ac + i + 1)
            if pi < pa:
                continue
            st_new = st + exp(pa - pi)
            if st_new == st:
                break
            st = st_new
    return max(0, pa - p0 - log(st))


# right tail only
def test_right(a, b, c, d):
    return exp(-mlnTest2r(a, a + b, a + c, a + b + c + d))


def mlnTest2r(a, ab, ac, abcd):
    if 0 > a or a > ab or a > ac or ab + ac > abcd + a:
        raise ValueError("invalid contingency table")
    if abcd > MAXN:
        raise OverflowError("the grand total of contingency table is too large")
    a_min = max(0, ab + ac - abcd)
    a_max = min(ab, ac)
    if a_min == a_max:
        return 0.0
    p0 = lgamma(ab + 1) + lgamma(ac + 1) + lgamma(abcd - ac + 1) + lgamma(abcd - ab + 1) - lgamma(abcd + 1)
    pa = lgamma(a + 1) + lgamma(ab - a + 1) + lgamma(ac - a + 1) + lgamma(abcd - ab - ac + a + 1)
    if ab * ac > a * abcd:
        sl = 0.0
        for i in range(a - 1, a_min - 1, -1):
            sl_new = sl + exp(
                pa - lgamma(i + 1) - lgamma(ab - i + 1) - lgamma(ac - i + 1) - lgamma(abcd - ab - ac + i + 1)
            )
            if sl_new == sl:
                break
            sl = sl_new
        return -log(1.0 - max(0, exp(p0 - pa) * sl))
    else:
        sr = 1.0
        for i in range(a + 1, a_max + 1):
            sr_new = sr + exp(
                pa - lgamma(i + 1) - lgamma(ab - i + 1) - lgamma(ac - i + 1) - lgamma(abcd - ab - ac + i + 1)
            )
            if sr_new == sr:
                break
            sr = sr_new
        return max(0, pa - p0 - log(sr))


# left tail only
def test_left(a, b, c, d):
    return exp(-mlnTest2l(a, a + b, a + c, a + b + c + d))


def mlnTest2l(a, ab, ac, abcd):
    if 0 > a or a > ab or a > ac or ab + ac > abcd + a:
        raise ValueError("invalid contingency table")
    if abcd > MAXN:
        raise OverflowError("the grand total of contingency table is too large")
    a_min = max(0, ab + ac - abcd)
    a_max = min(ab, ac)
    if a_min == a_max:
        return 0.0
    p0 = lgamma(ab + 1) + lgamma(ac + 1) + lgamma(abcd - ac + 1) + lgamma(abcd - ab + 1) - lgamma(abcd + 1)
    pa = lgamma(a + 1) + lgamma(ab - a + 1) + lgamma(ac - a + 1) + lgamma(abcd - ab - ac + a + 1)
    if ab * ac < a * abcd:
        sr = 0.0
        for i in range(a + 1, a_max + 1):
            sr_new = sr + exp(
                pa - lgamma(i + 1) - lgamma(ab - i + 1) - lgamma(ac - i + 1) - lgamma(abcd - ab - ac + i + 1)
            )
            if sr_new == sr:
                break
            sr = sr_new
        return -log(1.0 - max(0, exp(p0 - pa) * sr))
    else:
        sl = 1.0
        for i in range(a - 1, a_min - 1, -1):
            sl_new = sl + exp(
                pa - lgamma(i + 1) - lgamma(ab - i + 1) - lgamma(ac - i + 1) - lgamma(abcd - ab - ac + i + 1)
            )
            if sl_new == sl:
                break
            sl = sl_new
        return max(0, pa - p0 - log(sl))
