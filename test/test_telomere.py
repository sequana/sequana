import pytest

from sequana.telomere import Telomere, circular_shifts, factorize_sequences

from . import test_dir

TELOMERE_FA = f"{test_dir}/data/fasta/test_telomere.fasta"
MEASLES_FA = f"{test_dir}/data/fasta/measles.fa"


# ---------------------------------------------------------------------------
# Module-level utilities
# ---------------------------------------------------------------------------


def test_circular_shifts():
    shifts = circular_shifts("AACCCT")
    assert len(shifts) == 6
    assert shifts[0] == "AACCCT"
    assert shifts[1] == "ACCCTA"
    assert "TAACCC" in shifts


def test_factorize_sequences():
    seqs = ["AACCCT", "ACCCTA", "CCCTAA"]  # all rotations of each other
    groups = factorize_sequences(seqs)
    # All three should map to a single representative
    assert len(groups) == 1
    members = list(groups.values())[0]
    assert len(members) == 3


# ---------------------------------------------------------------------------
# Telomere class — construction and k-mer discovery
# ---------------------------------------------------------------------------


def test_telomere_init_no_file():
    t = Telomere()
    assert len(t.kmers) == 12  # 6 canonical + 6 additional
    assert t.peak_height == 20
    assert t.peak_width == 50


def test_find_representative_kmers_telomere():
    t = Telomere(TELOMERE_FA)
    candidates = t.find_representative_kmers(t.fasta.sequences[0])
    assert "TAACCC" in [x[0] for x in candidates]


def test_find_representative_kmers_no_telomere():
    t = Telomere(MEASLES_FA)
    candidates = t.find_representative_kmers(t.fasta.sequences[0])
    assert len(candidates) == 0


# ---------------------------------------------------------------------------
# Sliding-window counts
# ---------------------------------------------------------------------------


def test_sliding_window_counts():
    t = Telomere(TELOMERE_FA)
    # Use a short synthetic sequence so the test stays fast
    seq = "AACCCT" * 50 + "N" * 200 + "AGGGTT" * 50
    XX = t.get_sliding_kmer_count_five_to_three_prime(seq, W=50)
    YY = t.get_sliding_kmer_count_three_to_five_prime(seq, W=50)
    assert len(XX) == len(seq)
    assert len(YY) == len(seq)
    # Forward signal should peak at the start
    assert XX[:100].mean() > XX[200:].mean()
    # Reverse signal should peak at the end
    assert YY[-100:].mean() > YY[:200].mean()


def test_is_telomeric():
    t = Telomere(TELOMERE_FA)
    telomeric_seq = "AACCCT" * 100
    nontelo_seq = "ACGTACGTACGT" * 50
    assert t.is_telomeric(telomeric_seq) > t.is_telomeric(nontelo_seq)


# ---------------------------------------------------------------------------
# run() — default annotated style
# ---------------------------------------------------------------------------


def test_run_with_telomeres(tmpdir):
    t = Telomere(TELOMERE_FA)
    df = t.run(tag=None)
    assert "telomere" in df.columns
    assert "length" in df.columns
    assert "length_wo_telomere" in df.columns
    # The test FASTA has a contig with both telomeres
    assert (df["telomere"] == "complete").any()


def test_run_without_telomeres(tmpdir):
    t = Telomere(MEASLES_FA)
    df = t.run(tag=None)
    assert "telomere" in df.columns
    assert (df["telomere"] == "none").any()


def test_run_legacy_style():
    t = Telomere(MEASLES_FA)
    df = t.run(tag=None, plot_style="legacy")
    assert "telomere" in df.columns


def test_run_subset_chromosomes():
    t = Telomere(TELOMERE_FA)
    names = t.fasta.names[:1]
    df = t.run(tag=None, names=names)
    assert len(df) == 1


# ---------------------------------------------------------------------------
# plot_contig()
# ---------------------------------------------------------------------------


def test_plot_contig():
    import matplotlib

    matplotlib.use("Agg")
    import numpy as np

    t = Telomere(TELOMERE_FA)
    L = 500
    XX = np.zeros(L)
    YY = np.zeros(L)
    # Simulate a LHS peak on XX and an RHS peak on YY
    XX[:50] = 30
    YY[-50:] = 30

    fig = t.plot_contig(
        XX,
        YY,
        chrom="chr1",
        total_length=1_000_000,
        midpoint=0,
        lhs1=50,
        rhs1=0,
        lhs1_extend=50,
        rhs1_extend=0,
        lhs2=0,
        rhs2=50,
        lhs2_extend=0,
        rhs2_extend=50,
    )
    assert fig is not None
    import pylab

    pylab.close(fig)


def test_plot_contig_with_reversed():
    import matplotlib

    matplotlib.use("Agg")
    import numpy as np

    t = Telomere(TELOMERE_FA)
    L = 500
    XX = np.zeros(L)
    YY = np.zeros(L)
    XX[-50:] = 30  # reversed RHS on forward strand

    fig = t.plot_contig(
        XX,
        YY,
        chrom="chrRev",
        total_length=500_000,
        midpoint=0,
        lhs1=0,
        rhs1=50,
        lhs1_extend=0,
        rhs1_extend=50,
        lhs2=0,
        rhs2=0,
        lhs2_extend=0,
        rhs2_extend=0,
    )
    assert fig is not None
    import pylab

    pylab.close(fig)


def test_plot_contig_with_midpoint():
    """Covers the trimmed-sequence code path where midpoint < L."""
    import matplotlib

    matplotlib.use("Agg")
    import numpy as np

    t = Telomere(TELOMERE_FA)
    L = 600
    XX = np.zeros(L)
    YY = np.zeros(L)
    XX[:50] = 25
    YY[-50:] = 25

    fig = t.plot_contig(
        XX,
        YY,
        chrom="chrTrimmed",
        total_length=2_000_000,
        midpoint=300,
        lhs1=50,
        rhs1=0,
        lhs1_extend=50,
        rhs1_extend=0,
        lhs2=0,
        rhs2=50,
        lhs2_extend=0,
        rhs2_extend=50,
    )
    assert fig is not None
    import pylab

    pylab.close(fig)


# ---------------------------------------------------------------------------
# plot_summary()
# ---------------------------------------------------------------------------


def test_plot_summary():
    import matplotlib

    matplotlib.use("Agg")
    import pandas as pd

    t = Telomere(TELOMERE_FA)

    # Build a minimal DataFrame mimicking run() output
    df = pd.DataFrame(
        {
            "name": ["chr1", "chr2", "chr3", "chr4"],
            "length": [1_000_000, 800_000, 600_000, 400_000],
            "5to3_LHS_position": [5000, 0, 4000, 0],
            "5to3_LHS_length": [5000, 0, 4000, 0],
            "5to3_RHS_position": [0, 0, 4500, 0],
            "5to3_RHS_length": [0, 0, 4500, 0],
            "3to5_LHS_position": [0, 0, 0, 0],
            "3to5_LHS_length": [0, 0, 0, 0],
            "3to5_RHS_position": [4800, 5200, 0, 0],
            "3to5_RHS_length": [4800, 5200, 0, 0],
            "telomere": ["complete", "RHS_only", "none", "none"],
            "length_wo_telomere": [990200, 794800, 591500, 400_000],
        }
    )

    fig = t.plot_summary(df)
    assert fig is not None
    import pylab

    pylab.close(fig)
