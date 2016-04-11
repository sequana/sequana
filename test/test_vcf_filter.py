from sequana import vcf_filter
from . import data
pathdata = data.__path__[0]
from easydev import TempFile
import os
import filecmp


def test_vcf_filter():
    ft = TempFile()
    vcf_output_expected = pathdata + os.sep + "vcf_filter_output.vcf"
    vcf_record = vcf_filter.VCF(pathdata + os.sep + "vcf_filter_input.vcf")
    filter_dict = {"QUAL": 10000, "FREQ": 0.85, "INFO": {"DP": ">10",
                                                        "AO": ">200",
                                                        "SRP": "<100"}}
    vcf_record.filter_vcf(filter_dict, ft.name)
    compare_file = filecmp.cmp(ft.name, vcf_output_expected)
    assert compare_file == True
    ft.delete()


def test_constructor():
    vcf = vcf_filter.VCF('dummy')
    

def test_special_cases():
    from sequana import sequana_data
    data = sequana_data("test.vcf", "testing")
   

    # Test the different filter comparisons (< <=, >, >=)
    vcf = vcf_filter.VCF(data)
    filter_dict = {"QUAL": 10000, "FREQ": 0.5, 
             "INFO": {"DP": ">=10", "AO": ">200", "SRP": "<=100"}}
    for variant in vcf:
        vcf._filter_line(variant, filter_dict)

    # WRONG/MISSONG key in the filter INFO 
    vcf = vcf_filter.VCF(data)
    filter_dict = {"QUAL": 10000, "FREQ": 0.1, 
             "INFO": {"DP": ">=10", "AO": ">200", "DUMMY": "<=0"}}
    try:
        for variant in vcf:
            vcf._filter_line(variant, filter_dict)
        assert False
    except:
        assert True

    vcf = vcf_filter.VCF(data)
    filter_dict = {"QUAL": 10000, "FREQ": 0.5, 
             "INFO": {"DP": "<=10"}}
    for variant in vcf:
        vcf._filter_line(variant, filter_dict)
