from sequana import Homer

from . import test_dir


def test_homer():

    h = Homer(f"{test_dir}/data/annotatePeaks/example1.txt")
    h.pie_annotation()


def test_homer_empty():
    h = Homer(f"{test_dir}/data/annotatePeaks/example2_empty.txt")
    h.pie_annotation()
