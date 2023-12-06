from sequana.mh import MetropolisHasting


def test_mh():

    from sequana.mh import MetropolisHasting

    m = MetropolisHasting()
    m.Xtarget = [
        0.0,
        0.005,
        0.01,
        0.016,
        0.021,
        0.027,
        0.032,
        0.037,
        0.043,
        0.048,
        0.054,
        0.059,
        0.065,
        0.07,
        0.075,
        0.081,
        0.086,
        0.092,
        0.097,
        0.103,
    ]
    m.Ytarget = [83, 315, 611, 675, 1497, 5099, 7492, 2797, 842, 334, 117, 63, 33, 22, 11, 3, 3, 1, 0, 2]
    vec = m.simulate(100000)
    m.check(bins=100)
    m.diagnostics(bins=100)
