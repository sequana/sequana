import os

from sequana.cruciforms import Cruciforms

from . import test_dir


def test_crux(tmpdir):

    outname = tmpdir.join("test.bed")
    data = f"{test_dir}/data/fasta/measles.fa"

    c = Cruciforms(data)
    c.run()
    c.to_bed(outname)
