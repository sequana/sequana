import string

import numpy as np
import pandas as pd
import pylab
import pytest

from sequana.viz import corrplot


@pytest.fixture
def myinstance():
    letters = string.ascii_uppercase[0:10]
    df = pd.DataFrame(dict(((k, np.random.random(10) + ord(k) - 65) for k in letters)))
    klass = corrplot.Corrplot(df.corr())
    klass = corrplot.Corrplot(df)
    return klass


def test_correlation(myinstance):
    df1 = pd.DataFrame([[1, 2, 3, 4], [4, 5, 1, 2]])
    c1 = corrplot.Corrplot(df1)

    df2 = pd.DataFrame([[1, 2, 3, 4], [4, 5, 1, 2]]).corr()
    c2 = corrplot.Corrplot(df2)

    # in c1, the correlation is computed.
    assert (c1.df == c2.df).all().all() == True


def test_plot_square(myinstance):
    myinstance.plot(colorbar=False, method="square", shrink=0.9, rotation=45)


def test_binarise(myinstance):
    myinstance.plot(binarise_color=True, method="square")


def test_plot_text(myinstance):

    myinstance.plot()
    myinstance.plot(method="text", fontsize=8)


def test_cmap(myinstance):
    myinstance.plot(method="color", cmap="copper")

    # test fig parameter and cmap made of 3 colors
    fig = pylab.figure()
    myinstance.plot(fig=fig, method="color", cmap=("green", "orange", "red"))

    try:
        myinstance.plot(fig=fig, method="color", cmap=1)
        assert False
    except:
        assert True


def test_plot_color(myinstance):
    myinstance.plot(method="color")


def test_order(myinstance):
    myinstance.order(method="single", inplace=True)
    myinstance.plot(method="pie", shrink=0.8, order_metric="jaccard")


def test_figure(myinstance):
    myinstance.plot(fig=1)
    myinstance.plot(fig=None)
    fig = pylab.clf()
    myinstance.plot(fig=fig)

    ax = pylab.subplot(2, 1, 2)
    myinstance.plot(fig=fig, ax=ax)


def test_lower(myinstance):
    myinstance.plot(colorbar=False, method="circle", shrink=0.9, lower="circle")


def test_upper(myinstance):
    myinstance.plot(colorbar=False, method="circle", shrink=0.9, upper="circle")


def test_both(myinstance):
    myinstance.plot(colorbar=False, method="circle", shrink=0.9, upper="circle", lower="ellipse")


def test_wrong_method(myinstance):
    try:
        myinstance.plot(colorbar=False, method="dummy")
        assert False
    except ValueError:
        assert True
