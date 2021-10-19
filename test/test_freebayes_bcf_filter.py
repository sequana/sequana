import filecmp

from sequana.freebayes_bcf_filter import BCF_freebayes
from . import test_dir
sharedir = f"{test_dir}/data/vcf/"

def test_bcf_filter(tmpdir):
    path = tmpdir.mkdir("temp")
    vcf_output_expected = f"{sharedir}/JB409847.filter.vcf"


    bcf = BCF_freebayes(f"{sharedir}/JB409847.bcf")
    filter_dict = {'freebayes_score': 200, 'frequency': 0.85, 'min_depth': 10,
                   'forward_depth': 3, 'reverse_depth': 3, 'strand_ratio': 0.3}
    filter_bcf = bcf.filter_bcf(filter_dict)

    with open(path+"/test.vcf", "w") as fp:
        filter_bcf.to_vcf(fp.name)
        compare_file = filecmp.cmp(fp.name, vcf_output_expected)
        assert compare_file


def test_constructor():
    try:
        BCF_freebayes('dummy')
        assert False
    except OSError:
        assert True


def test_to_csv(tmpdir):
    path = tmpdir.mkdir("temp")

    filter_dict = {'freebayes_score': 200, 'frequency': 0.85, 'min_depth': 10,
                   'forward_depth': 3, 'reverse_depth': 3, 'strand_ratio': 0.3}
    bcf = BCF_freebayes(f"{sharedir}/JB409847.bcf")
    filter_bcf = bcf.filter_bcf(filter_dict)
    with open(path + "/test.csv", "w") as ft:
        filter_bcf.to_csv(ft.name)
