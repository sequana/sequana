from sequana import annotation

from . import test_dir


def test_aragorn_parser():

    arg = annotation.Aragorn()
    df = arg.parse_aragorn_output(f"{test_dir}/data/aragorn.txt")
    assert len(df) == 90
