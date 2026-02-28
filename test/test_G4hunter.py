import os

from easydev import md5 as md5sum
from numpy import mean

from sequana import FastA
from sequana.G4hunter import G4Hunter, G4HunterReader

from . import test_dir


def test_g4hunter_ecoli(tmpdir):
    outdir = tmpdir.mkdir("temp")

    out1 = outdir.join("ecoli_MG1655-W20-S1.txt")
    out2 = outdir.join("ecoli_MG1655-Merged.txt")

    filename = f"{test_dir}/data/fasta/ecoli_MG1655.fa"
    ff = G4Hunter(filename, 20, 1)
    ff.run(outdir)

    assert md5sum(out2) == "a135ccda22a3e398d3151d856a1dd6db"
    assert md5sum(out1) == "b59779479108de459a38033046c14afd"


def test_base_score():
    filename = f"{test_dir}/data/fasta/ecoli_MG1655.fa"
    ff = G4Hunter(filename, 20, 1)
    assert mean(ff.base_score("CCC")) == -3.0
    assert mean(ff.base_score("GGGGGG")) == 4.0
    assert mean(ff.base_score("ATCCCAAGGGAA")) == 0.0
    assert mean(ff.base_score("ATGGATGGATGATGAT")) == 0.625


def test_load_merged_data(tmpdir):

    # Create a temporary file with test data
    test_file = tmpdir.join("test_merged.txt")
    test_file.write(
        ">seq1\n"
        "Start\tEnd\tSequence\tLength\tScore\n"
        "0\t20\tACGTACGTACGTACGTACGT\t20\t1.5\n"
        "30\t50\tGGGGCCCCGGGGCCCC\t20\t2.0\n"
        ">seq2\n"
        "Start\tEnd\tSequence\tLength\tScore\tNBR\n"
        "10\t30\tATATATATATATATAT\t20\t0.8\t2\n"
        "40\t60\tCCCCGGGGCCCCGGGG\t20\t1.2\t3\n"
    )

    reader = G4HunterReader()
    reader.load_merged_data(str(test_file))

    assert "seq1" in reader.data_merged
    assert "seq2" in reader.data_merged
    assert len(reader.data_merged["seq1"]) == 2
    assert len(reader.data_merged["seq2"]) == 2

    # Verify seq1 data
    assert reader.data_merged["seq1"].iloc[0]["Start"] == 0
    assert reader.data_merged["seq1"].iloc[0]["End"] == 20
    assert reader.data_merged["seq1"].iloc[0]["Score"] == 1.5
    assert reader.data_merged["seq1"].iloc[0]["NBR"] is None

    # Verify seq2 data with NBR column
    assert reader.data_merged["seq2"].iloc[0]["Start"] == 10
    assert reader.data_merged["seq2"].iloc[0]["NBR"] == 2
    assert reader.data_merged["seq2"].iloc[1]["NBR"] == 3


def test_load_merged_data_empty_file(tmpdir):

    test_file = tmpdir.join("empty.txt")
    test_file.write("")

    reader = G4HunterReader()
    reader.load_merged_data(str(test_file))

    assert len(reader.data_merged) == 0


def test_load_merged_data_single_sequence(tmpdir):

    test_file = tmpdir.join("single.txt")
    test_file.write(">seq_single\n" "Start\tEnd\tSequence\tLength\tScore\n" "5\t25\tGGGGGCCCCCGGGGGCCCCC\t20\t2.5\n")

    reader = G4HunterReader()
    reader.load_merged_data(str(test_file))

    assert "seq_single" in reader.data_merged
    assert len(reader.data_merged["seq_single"]) == 1
    assert reader.data_merged["seq_single"].iloc[0]["Score"] == 2.5
