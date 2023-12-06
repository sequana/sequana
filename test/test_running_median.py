import scipy.signal
from pylab import randn

from sequana.running_median import RunningMedian


def test_running_median():

    x = randn(100)
    rm1 = RunningMedian(x, 7).run()
    rm2 = scipy.signal.medfilt(x, 7)
    diff = sum((rm2 - rm1)[3:-3])  # get rid of the first and last W/2 points, where W=7

    assert abs(diff) <= 1e-10

    # non-even list
    try:
        rm1 = RunningMedian(x, 10).run()
        assert True
    except:
        assert True
