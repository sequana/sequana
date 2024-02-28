import filecmp

import pytest

from sequana.freebayes import VCF_freebayes

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


def test_to_csv(tmpdir):
    path = tmpdir.mkdir("temp")

    filter_dict = {
        "freebayes_score": 200,
        "frequency": 0.85,
        "min_depth": 20,
        "forward_depth": 3,
        "reverse_depth": 3,
        "strand_ratio": 0.3,
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
