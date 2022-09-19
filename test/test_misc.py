from sequana.misc import *
from math import pi, sqrt


def test_normpdf():
    sigma = 2
    assert sqrt(2*pi)*sigma * normpdf(0, 0, sigma) == 1

def test_textwrap():
    res = textwrap("test1test2", width=5, indent=0)
    assert res.split("\n")[0] == "test1"
    assert res.split("\n")[1] == "test2"
    res = textwrap("test1test2", width=5, indent=4)


def test_findpos():
    assert list(findpos("AACCTTGGAACCGG", "GG"))


def test_wget():
    from easydev import TempFile
    with TempFile() as fh:
        wget("https://github.com/sequana/sequana/raw/main/README.rst", fh.name)

