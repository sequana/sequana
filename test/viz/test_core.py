import random

from sequana.viz import Hist2D, Imshow


def test1():
    # test input a list of 2 lists
    Imshow([[1, 2], [3, 4]])
    from pylab import randn

    X = randn(10000)
    Y = randn(10000)
    h = Hist2D(X, Y)
    h.plot(bins=100, contour=True)

    import pandas as pd

    df = pd.DataFrame({"A": [1, 2], "B": [3, 4]})
    h = Hist2D(df)

    h = Hist2D(x=[1, 2], y=[3, 4])
    h = Hist2D(x=[1, 2])
