from sequana import Homer

from . import test_dir

def test_homer():

    h = Homer(f"{test_dir}/data/annotatePeaks/example1.txt")
    h.pie_annotation()
