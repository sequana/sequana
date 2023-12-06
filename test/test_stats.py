from sequana.stats import evenness, moving_average


def test_ma():
    ma = moving_average([1, 1, 1, 1, 3, 3, 3, 3], 4)
    assert list(ma) == [1.0, 1.5, 2.0, 2.5, 3.0]


def test_evenness():
    assert evenness([1, 1, 1, 1, 4, 4, 4, 4]) == 0.75
    assert evenness([1, 1, 1, 1]) == 1


def test_runmean():
    from sequana.stats import runmean

    ma = runmean([1, 1, 1, 1, 3, 3, 3, 3], 4)
    assert list(ma) == [1.0, 1.5, 2.0, 2.5, 3.0]


def test_N50():
    from sequana.stats import N50

    assert 1002 == N50([1000, 1001, 1002, 1003])


def test_L50():
    from sequana.stats import L50

    assert 2 == L50([1000, 1001, 1002, 1003])
