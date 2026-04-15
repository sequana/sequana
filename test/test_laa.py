import os
import shutil

import matplotlib

matplotlib.use("Agg")

from sequana.laa import LAA, Consensus, LAA_Assembly

from . import test_dir

CSV_BC24 = os.path.join(test_dir, "data", "csv", "test_laa_amplicon_analysis_summary_bc24.csv")
CSV_BC25 = os.path.join(test_dir, "data", "csv", "test_laa_amplicon_analysis_summary_bc25.csv")
BAM = os.path.join(test_dir, "data", "bam", "test_measles.bam")


def _make_bases_file(path, npos=50, header_lines=3):
    with open(path, "w") as fout:
        fout.write("# header chrom=ref1\n")
        fout.write("# col Pos A C G T N DEL INS\n")
        for _ in range(header_lines - 2):
            fout.write("# extra header\n")
        bases = ["A", "C", "G", "T"]
        for pos in range(npos):
            counts = [0] * 7
            counts["ACGT".index(bases[pos % 4])] = 100
            row = "\t".join(str(v) for v in [pos] + counts)
            fout.write(row + "\n")


def test_laa_class(tmp_path):
    for name, src in [("bc24", CSV_BC24), ("bc25", CSV_BC25)]:
        d = tmp_path / name
        d.mkdir()
        shutil.copy(src, d / "amplicon_analysis_summary.csv")

    cwd = os.getcwd()
    try:
        os.chdir(tmp_path)
        laa = LAA(where="bc*")
        assert len(laa.filenames) == 2
        assert len(laa.data) == 2
        assert laa.numbers == [len(d) for d in laa.data]
        laa.hist_amplicon()
        laa.plot_max_length_amplicon_per_barcode()
        laa.plot_max_length_amplicon_per_barcode(sample_names=["a", "b"])
    finally:
        os.chdir(cwd)


def test_consensus_get_bases(tmp_path):
    bases = tmp_path / "bases.txt"
    _make_bases_file(str(bases), npos=20)
    c = Consensus(str(bases))
    df = c.get_bases()
    assert list(df.columns) == ["A", "C", "G", "T", "N", "DEL", "INS"]
    assert len(df) == 20


def test_consensus_get_bases_low_depth(tmp_path):
    bases = tmp_path / "bases.txt"
    with open(bases, "w") as fout:
        fout.write("# header chrom=ref1\n")
        fout.write("# col\n")
        fout.write("# col\n")
        fout.write("0\t1\t0\t0\t0\t0\t0\t0\n")
        fout.write("1\t100\t0\t0\t0\t0\t0\t0\n")
    c = Consensus(str(bases))
    df = c.get_bases()
    assert df.loc[0, "N"] == 10000
    assert df.loc[0, "A"] == 0
    assert df.loc[1, "A"] == 100


def test_consensus_run_and_save(tmp_path):
    bases = tmp_path / "bases.txt"
    _make_bases_file(str(bases), npos=20)
    c = Consensus(str(bases))
    seq = c.run()
    assert len(seq) == 20
    assert "".join(seq[:4]) == "ACGT"

    out = tmp_path / "consensus.fa"
    c.save_consensus(str(out), "sample1")
    content = out.read_text()
    assert content.startswith(">sample1\n")
    assert "ACGT" in content


def test_consensus_get_population(tmp_path):
    bases = tmp_path / "bases.txt"
    _make_bases_file(str(bases), npos=10)
    c = Consensus(str(bases))
    sel = c.get_population(threshold=0.1, Npop=2)
    assert len(sel) == 0


def test_consensus_identify_deletions_no_freebayes(tmp_path):
    bases = tmp_path / "bases.txt"
    _make_bases_file(str(bases), npos=5)
    c = Consensus(str(bases))
    assert c.identify_deletions() == []


def test_consensus_split_according_to_reference_name(tmp_path):
    bases = tmp_path / "input.txt"
    with open(bases, "w") as fout:
        fout.write("# global header line 1\n")
        fout.write("# global header line 2\n")
        fout.write("# section chrom=chrA blah\n")
        fout.write("0\t100\t0\t0\t0\t0\t0\t0\n")
        fout.write("1\t100\t0\t0\t0\t0\t0\t0\n")
        fout.write("# section chrom=chrB blah\n")
        fout.write("0\t0\t100\t0\t0\t0\t0\t0\n")
    c = Consensus(str(bases))
    cwd = os.getcwd()
    try:
        os.chdir(tmp_path)
        c.split_according_to_reference_name(prefix="cons")
        assert os.path.exists("cons_chrA.txt")
        assert os.path.exists("cons_chrB.txt")
    finally:
        os.chdir(cwd)


def test_laa_assembly(tmp_path):
    asm = LAA_Assembly(BAM)
    seq = asm.build_reference()
    assert isinstance(seq, str)
    assert len(seq) > 0

    out = tmp_path / "out.fa"
    asm.save_fasta(str(out), sequence="ACGT")
    assert out.read_text().startswith(">test\nACGT")

    out2 = tmp_path / "out2.fa"
    asm.save_fasta(str(out2))
    assert out2.read_text().startswith(">test\n")
