import os

from sequana import FastA
from sequana.G4hunter import G4Hunter
from easydev import md5 as md5sum
from numpy import mean

from . import test_dir


def test_g4hunter(tmpdir):
    outdir = tmpdir.mkdir("temp")


    out1 = outdir.join('ecoli_MG1655-W20-S1.txt')
    out2 = outdir.join('ecoli_MG1655-Merged.txt')

    filename = f"{test_dir}/data/fasta/ecoli_MG1655.fa"
    ff = G4Hunter(filename, 20, 1)
    ff.run(outdir)

    assert md5sum(out2) == "a135ccda22a3e398d3151d856a1dd6db"
    assert md5sum(out1) == "b59779479108de459a38033046c14afd"

    assert mean(ff.base_score_python("CCC")) == -3.0
    assert mean(ff.base_score_python("GGGGGG")) == 4.0
    assert mean(ff.base_score_python("ATCCCAAGGGAA")) == 0.0
    assert mean(ff.base_score_python("ATGGATGGATGATGAT")) == 0.625

