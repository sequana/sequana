from sequana import GFF3, FastA
from sequana.kozak import Kozak

from . import test_dir


def test_kozak():

    fastafile = f"{test_dir}/data/fasta/ecoli_MG1655.fa"
    gff = f"{test_dir}/data/gff/ecoli_MG1655.gff"

    k = Kozak(fastafile, gff, "gene", "ID")
    k.set_context(left_kozak=15, right_kozak=12, keep_ATG_only=True, include_start_codon=True)

    df = k.get_data()

    assert len(df["kozak_left"].iloc[0]) == 15
    assert list(df["start_codon"].unique()) == ["ATG"]
    assert k.metrics["ATG_ratio"] != 1  # original data has non-ATG start codon

    # check caching is functional (correctly reset)
    k.set_context(left_kozak=12, right_kozak=12, keep_ATG_only=False, include_start_codon=True)
    df = k.get_data()
    assert len(df["kozak_left"].iloc[0]) == 12
    assert k.metrics["ATG_ratio"] != 1
    assert list(df["start_codon"].unique()) != ["ATG"]

    dd = k.plot_logo(df)
    k.plot_logo_purine_pyrimidine(df)
    k.get_gc_per_chromosome()

    assert k._get_logo_data(df).equals(k.plot_logo())

    motif = k._get_logo_data(df)
    assert k._add_purine_pyrimidine(motif).equals(k.plot_logo_purine_pyrimidine())

    k.get_entropy(motif)
