import os

from easydev import TempFile
from pylab import close as _plt_close

from sequana import FastA
from sequana.fasta import is_fasta

from . import test_dir


def test_lazy_sequences_basic():
    filename = f"{test_dir}/data/test_fasta.fasta"
    ff = FastA(filename)

    # behaves like list
    assert len(ff.sequences) == len(ff.names)

    # indexing
    first = ff.sequences[0]
    second = ff.sequences[1]
    assert isinstance(first, str)
    assert isinstance(second, str)
    assert len(first) > 0

    # iteration
    seqs = list(ff.sequences)
    assert len(seqs) == len(ff.names)


def test_lazy_sequences_slice():
    filename = f"{test_dir}/data/test_fasta.fasta"
    ff = FastA(filename)

    subset = ff.sequences[2:5]
    assert isinstance(subset, list)
    assert len(subset) == 3
    for s in subset:
        assert isinstance(s, str)


def test_lazy_sequences_repeat_access():
    filename = f"{test_dir}/data/test_fasta.fasta"
    ff = FastA(filename)

    # access same index twice (cache test)
    a = ff.sequences[3]
    b = ff.sequences[3]
    assert a == b


def test_lazy_sequences_iter_partial():
    filename = f"{test_dir}/data/test_fasta.fasta"
    ff = FastA(filename)

    # partial iteration should work without exhausting object
    it = iter(ff.sequences)
    first = next(it)
    second = next(it)

    assert isinstance(first, str)
    assert isinstance(second, str)


def test_lazy_cache_clear():
    filename = f"{test_dir}/data/test_fasta.fasta"
    ff = FastA(filename)

    _ = ff.sequences[0]
    _ = ff.sequences[1]

    # clear cache should not break behavior
    ff.sequences.clear_cache()

    again = ff.sequences[0]
    assert isinstance(again, str)


def test_format_contigs_denovo():
    # test with a custom fasta
    filename = f"{test_dir}/data/test_fasta.fasta"
    contigs = FastA(filename)
    with TempFile(suffix=".fasta") as fh:
        contigs.format_contigs_denovo(fh.name)
    contigs.names
    contigs.lengths
    contigs.comments
    contigs.GC_content()
    import pytest

    with pytest.warns(DeprecationWarning):
        contigs.GC_content_sequence(contigs._fasta[contigs.names[0]])

    gc = contigs.GC_content_per_sequence()
    assert sorted(gc.keys()) == sorted(contigs.names)
    assert all(0 <= v <= 100 for v in gc.values())

    excluded = contigs.names[0]
    gc2 = contigs.GC_content_per_sequence(exclude=[excluded])
    assert excluded not in gc2
    assert len(gc2) == len(gc) - 1

    contigs.plot_GC_content()
    contigs.plot_GC_content(kind="bar", max_contigs=1, sorting="alphanum", hline=False)
    contigs.plot_GC_content(exclude=[excluded])
    _plt_close("all")

    contigs.summary()


def test_fasta_filtering():
    filename = f"{test_dir}/data/test_fasta_filtering.fasta"
    ff = FastA(filename)
    with TempFile(suffix=".fasta") as fh:
        ff.to_fasta(fh.name)
        ff.save_ctg_to_fasta("A", fh.name)
        ff.save_ctg_to_fasta("A", fh.name, max_length=1)

    with TempFile(suffix=".fasta") as fh:
        ff.filter(fh.name, names_to_exclude=["A", "B"])
        reader = FastA(fh.name)
        assert set(reader.names) == set(["C", "D"])

    ff = FastA(filename)
    with TempFile(suffix=".fasta") as fh:
        ff.filter(
            fh.name,
            names_to_keep=[
                "A",
            ],
        )
        reader = FastA(fh.name)
        assert set(reader.names) == set(["A"])


def test_sequences_not_list():
    filename = f"{test_dir}/data/test_fasta.fasta"
    ff = FastA(filename)

    assert not isinstance(ff.sequences, list)


def test_others():
    filename = f"{test_dir}/data/test_fasta.fasta"
    ff = FastA(filename)
    assert len(ff) == 16
    assert len(ff.comments) == 16
    assert len(ff.names) == 16
    assert len(ff.sequences) == 16
    assert is_fasta(filename) == True
    ff.get_lengths_as_dict()
    with TempFile(suffix=".fasta") as fh:
        ff.select_random_reads(4, output_filename=fh.name)
        ff.select_random_reads([1, 2, 3], output_filename=fh.name)
        ff.select_random_reads({1, 2, 3}, output_filename=fh.name)
        ff.select_random_reads(100000, output_filename=fh.name)
    assert ff.get_stats()["N"] == 16
    assert ff.get_stats()["mean_length"] > 454
    with TempFile(suffix=".fasta") as fh:
        ff.reverse_and_save(fh.name)
        ff.to_fasta(fh.name)
        ff.to_fasta(fh.name)
        ff.to_igv_chrom_size(fh.name)

    with TempFile(suffix=".fasta") as fh:
        ff.save_collapsed_fasta(fh.name, "main")
        ff.save_collapsed_fasta(fh.name, "main", comment="null")

    ff.get_cumulative_sum("mixed")
    ff.get_cumulative_sum("alphanum")
    ff.find_gaps()


def test_explode(tmpdir):
    path = tmpdir.mkdir("temp")

    filename = f"{test_dir}/data/test_fasta.fasta"
    ff = FastA(filename)
    with TempFile(suffix=".fasta") as fh:
        ff.explode(outdir=path)


def test_fasta_comparison(tmp_path):
    from sequana.fasta import FastAComparison

    seq1 = "ACGTACGTTTGACCATGCA" * 5
    seq2 = "GGGTTTCCCAAA" * 4

    f1 = tmp_path / "f1.fa"
    f1.write_text(f">1\n{seq1}\n>2\n{seq2}\n")
    f2 = tmp_path / "f2.fa"
    f2.write_text(f">NC_0001.1\n{seq1}\n>NC_0002.1\n{seq2}\n")

    c = FastAComparison(f1, f2)
    assert c.identical is True
    assert c.orphans1 == []
    assert c.orphans2 == []
    assert c.matches == [(["1"], ["NC_0001.1"], len(seq1), False), (["2"], ["NC_0002.1"], len(seq2), False)]

    # a missing checksum must not be added to the internal dictionaries
    assert len(c.data1) == 2
    assert len(c.data2) == 2


def test_fasta_comparison_orphans_and_duplicates(tmp_path):
    from sequana.fasta import FastAComparison

    seq1 = "ACGTACGTTTGACCATGCA" * 5
    seq2 = "GGGTTTCCCAAA" * 4

    f1 = tmp_path / "f1.fa"
    f1.write_text(f">1\n{seq1}\n>1_bis\n{seq1}\n")
    f2 = tmp_path / "f2.fa"
    f2.write_text(f">NC_0001.1\n{seq1}\n>NC_0002.1\n{seq2}\n")

    c = FastAComparison(f1, f2)
    assert c.identical is False
    assert c.orphans1 == []
    assert c.orphans2 == [("NC_0002.1", len(seq2))]
    assert c.duplicates1 == [["1", "1_bis"]]
    assert c.duplicates2 == []
    names1, names2, length, rc = c.matches[0]
    assert names1 == ["1", "1_bis"]
    assert names2 == ["NC_0001.1"]
    assert rc is False


def test_fasta_comparison_case_gaps_and_rc(tmp_path):
    from sequana.fasta import FastAComparison
    from sequana.tools import reverse_complement

    seq = "ACGTACGTTTGACCATGCA" * 5

    f1 = tmp_path / "f1.fa"
    f1.write_text(f">1\n{seq}\n")
    f2 = tmp_path / "f2.fa"
    f2.write_text(f">NC_0001.1\n{seq.lower()}\n")
    assert FastAComparison(f1, f2).identical is True
    assert FastAComparison(f1, f2, ignore_case=False).identical is False

    f3 = tmp_path / "f3.fa"
    f3.write_text(f">NC_0001.1\n{seq[:10]}--{seq[10:]}\n")
    assert FastAComparison(f1, f3).identical is False
    assert FastAComparison(f1, f3, ignore_gaps=True).identical is True

    f4 = tmp_path / "f4.fa"
    f4.write_text(f">NC_0001.1\n{reverse_complement(seq)}\n")
    assert FastAComparison(f1, f4).identical is False
    c = FastAComparison(f1, f4, rc_aware=True)
    assert c.identical is True
    assert c.matches[0][3] is True
