import filecmp

import pytest

from sequana.freebayes_vcf_filter import VCF_freebayes

# just to import this alias
from sequana.variants import VariantFile

from . import test_dir

sharedir = f"{test_dir}/data/vcf/"


def test_vcf_filter(tmpdir):
    path = tmpdir.mkdir("temp")

    vcf_output_expected = f"{sharedir}/JB409847.expected.vcf"
    v = VCF_freebayes(f"{sharedir}/JB409847.vcf")
    filter_dict = {
        "freebayes_score": 200,
        "frequency": 0.85,
        "min_depth": 10,
        "forward_depth": 3,
        "reverse_depth": 3,
        "strand_ratio": 0.3,
        "keep_polymorphic": True,
    }
    filter_v = v.filter_vcf(filter_dict)
    assert len(filter_v.variants) == 24
    with open(path + "/test.vcf", "w") as ft:
        filter_v.to_vcf(ft.name)
        compare_file = filecmp.cmp(ft.name, vcf_output_expected)
        assert compare_file

    v.barplot()
    v.manhattan_plot("JB409847")


def test_constructor():
    with pytest.raises(OSError):
        VCF_freebayes("dummy")


def test_empty_vcf(tmpdir):
    # header only, no variant: must not raise StopIteration
    v = VCF_freebayes(f"{sharedir}/empty.vcf")
    assert len(v) == 0
    assert v._snpeff is False
    assert v.df.shape == (0, 0)

    filter_dict = {
        "freebayes_score": 20,
        "frequency": 0.1,
        "min_depth": 10,
        "forward_depth": 3,
        "reverse_depth": 3,
        "strand_ratio": 0.2,
        "keep_polymorphic": True,
    }
    filter_v = v.filter_vcf(filter_dict)
    assert len(filter_v.variants) == 0


def test_to_csv(tmpdir):
    path = tmpdir.mkdir("temp")

    filter_dict = {
        "freebayes_score": 200,
        "frequency": 0.85,
        "min_depth": 20,
        "forward_depth": 3,
        "reverse_depth": 3,
        "strand_ratio": 0.3,
        "keep_polymorphic": True,
    }
    v = VCF_freebayes(f"{sharedir}/JB409847.expected.vcf")
    filter_v = v.filter_vcf(filter_dict)
    assert len(filter_v.variants) == 3

    with open(path + "/test.csv", "w") as ft:
        filter_v.to_csv(ft.name)


def test_variant():
    v = VCF_freebayes(f"{sharedir}/JB409847.vcf")
    variants = v.variants
    assert len(variants) == 64
    print(variants[0])

    v = VCF_freebayes(f"{sharedir}/test_vcf_snpeff.vcf")
    variants = v.variants
    assert len(variants) == 775


def test_get_variant_type():
    v = VariantFile(f"{sharedir}/JB409847.vcf")
    variant_types = v.get_variant_type()
    assert variant_types["snp"] == 62
    assert variant_types["complex"] == 2


def test_barplot_and_pieplot():
    v = VariantFile(f"{sharedir}/JB409847.vcf")
    v.barplot()
    v.pieplot()


def test_manhattan_plot():
    v = VariantFile(f"{sharedir}/JB409847.vcf")
    v.manhattan_plot()
    v.manhattan_plot(chrom_name="JB409847")


def test_joint_calling_multi_sample_and_snpeff():
    """Joint calling VCF: multiple samples (is_joint branch) + snpEff EFF
    annotation, exercising the multi-sample frequency/GL/effect-parsing code
    paths not covered by the single-sample fixtures above."""
    v = VariantFile(f"{sharedir}/joint_calling.vcf")
    assert v.is_joint is True
    assert len(v.samples) == 5

    df = v.df
    assert len(df) == 7
    assert "effect_type" in df.columns
    assert "gene_name" in df.columns
    # one info_N column per sample
    for i in range(5):
        assert f"info_{i}" in df.columns


def test_filtered_variant_file_to_csv_with_info_field(tmpdir):
    path = tmpdir.mkdir("temp")
    v = VariantFile(f"{sharedir}/joint_calling.vcf")
    filtered = v.filter_vcf(
        {
            "freebayes_score": 0,
            "frequency": 0,
            "min_depth": 0,
            "forward_depth": 0,
            "reverse_depth": 0,
            "strand_ratio": 0,
            "keep_polymorphic": True,
        }
    )
    out = str(path.join("with_info.csv"))
    filtered.to_csv(out, info_field=True)
    with open(out) as f:
        content = f.read()
    assert "info_0" in content


def test_apply_variants_applies_snp(tmpdir):
    """apply_variants() must correctly substitute a known SNP: JB409847 VCF
    has a C->T substitution at position 2221 (1-based)."""
    from sequana.variants import apply_variants

    out_fasta = str(tmpdir.join("consensus.fasta"))
    apply_variants(f"{sharedir}/JB409847.fasta", f"{sharedir}/JB409847.vcf", out_fasta)

    import pysam

    original = pysam.FastaFile(f"{sharedir}/JB409847.fasta").fetch("JB409847")
    assert original[2220] == "C"

    consensus = pysam.FastaFile(out_fasta).fetch("JB409847")
    assert consensus[2220] == "T"
    # sequence length is preserved for a pure SNP substitution
    assert len(consensus) == len(original)


def test_variant_file_iterator_protocol():
    """VariantFile implements __iter__/__len__/__next__ so it can be looped
    over directly (e.g. `for v in VariantFile(...)`), independently of the
    `.variants` property."""
    v = VariantFile(f"{sharedir}/JB409847.vcf")
    assert len(v) == 64

    collected = list(v)
    assert len(collected) == 64

    # iterator resets after exhaustion (StopIteration resets the index)
    collected_again = list(v)
    assert len(collected_again) == 64


def test_hist_score():
    v = VariantFile(f"{sharedir}/JB409847.vcf")
    v.hist_score()
    v.hist_score(min_score=100)
