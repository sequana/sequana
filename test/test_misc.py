from math import pi, sqrt

from sequana.misc import *


def test_normpdf():
    sigma = 2
    assert sqrt(2 * pi) * sigma * normpdf(0, 0, sigma) == 1


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


def test_multiple_downloads(tmpdir): 
    file1 = tmpdir.join("file1.txt")
    file2 = tmpdir.join("file2.txt")
    data = [
        ("https://raw.githubusercontent.com/sequana/sequana_pipetools/main/README.rst", file1, 0),
        ("https://raw.githubusercontent.com/sequana/sequana_pipetools/main/requirements.txt", file2, 1),
    ]
    multiple_downloads(data)
    download(data[0][0], data[0][1])

