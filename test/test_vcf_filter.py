from sequana.vcf_filter import VCF
from sequana import sequana_data
from easydev import TempFile




def test_vcf_filter_freebayes():

    data = sequana_data("test.vcf")
    v = VCF(data)
    v.vcf.hist_qual()


def test_vcf_filter():

    data = sequana_data("test_vcf_mpileup_4dot1.vcf")
    v = VCF(data)
    assert v.vcf.version == "4.1"


    # Test the INFO set to {}
    v.vcf.filter_dict['INFO'] = {}
    v.vcf.filter_dict['QUAL'] = 50
    with TempFile() as fh:
        res = v.vcf.filter_vcf(fh.name)
    assert res == {'N': 573, 'filtered': 308, 'unfiltered': 265}


    # Test the & filter
    v.vcf.filter_dict["QUAL"] = 0
    v.vcf.filter_dict["INFO"] = {}
    v.vcf.rewind();
    v.vcf.filter_dict['INFO']['DP'] = ">40&<=80"
    with TempFile() as fh:
        res = v.vcf.filter_vcf(fh.name)
    assert res == {'N': 573, 'filtered': 10, 'unfiltered': 563}


    # check the | filter
    v.vcf.filter_dict["QUAL"] = 0
    v.vcf.filter_dict["INFO"] = {}
    v.vcf.rewind();
    v.vcf.filter_dict['INFO']['DP'] = "<40|>=80"
    with TempFile() as fh:
        res = v.vcf.filter_vcf(fh.name)
    assert res ==  {'N': 573, 'filtered': 562, 'unfiltered': 11}


    # check the array filter
    v.vcf.filter_dict["QUAL"] = 0
    v.vcf.filter_dict["INFO"] = {}
    v.vcf.rewind();
    v.vcf.filter_dict['INFO']['DP4[0]'] = "<2"
    with TempFile() as fh:
        res = v.vcf.filter_vcf(fh.name)
    assert res ==  {'N': 573, 'filtered': 173, 'unfiltered': 400}

    # Check the array access to one item
    v.vcf.filter_dict["QUAL"] = 0
    v.vcf.filter_dict["INFO"] = {}
    v.vcf.rewind();
    v.vcf.filter_dict['INFO']['DP4[2]'] = "<2"
    with TempFile() as fh:
        res = v.vcf.filter_vcf(fh.name)
    assert res ==  {'N': 573, 'filtered': 199, 'unfiltered': 374}

    # Check the sum of elements in an array
    v.vcf.filter_dict["QUAL"] = 0
    v.vcf.filter_dict["INFO"] = {}
    v.vcf.rewind();
    v.vcf.filter_dict['INFO']['sum(DP4[2], DP4[3])'] = "<4"
    with TempFile() as fh:
        res = v.vcf.filter_vcf(fh.name)
    assert res ==  {'N': 573, 'filtered': 7, 'unfiltered': 566}


def test_af1():
    data = sequana_data("test_vcf_mpileup_4dot1.vcf")
    v = VCF(data)
    variant = next(v.vcf)

    variant.INFO['AF1'] = 1
    assert v.vcf.is_valid_af1(variant) is True
    variant.INFO['AF1'] = 0.5
    assert v.vcf.is_valid_af1(variant) is False

    # polymorphic case
    variant = next(v.vcf)

    variant.INFO['AF1'] = 1
    assert v.vcf.is_valid_af1(variant) is False
    #variant.INFO['AF1'] = 0.5
    #assert v.vcf.is_valid_af1(variant) is True

    v = VCF(data)
    v.vcf.apply_af1_filter = True
    with TempFile() as fh:
        res = v.vcf.filter_vcf(fh.name)
    assert res ==  {'N': 573, 'filtered': 391, 'unfiltered': 182}


def test_indel():
    data = sequana_data("test_vcf_mpileup_4dot1.vcf")
    v = VCF(data)
    variant = next(v.vcf)
    assert v.vcf.is_indel(variant) is False
    variant = next(v.vcf)
    assert v.vcf.is_indel(variant) is False
    variant = next(v.vcf)
    assert v.vcf.is_indel(variant) is True


def test_vcf_filter_dp4():

    data = sequana_data("test_vcf_mpileup_4dot1.vcf")
    v = VCF(data)
    variant = next(v.vcf)

    def validate_variant_alternatate(variant):
        # variant.ALT must be different from "." for this test
        assert str(variant.ALT[0]).strip() != "."

        # test minimum depth of alternate must be >= 4
        variant.INFO['DP4'] = [0,0,2,2]
        assert v.vcf.is_valid_dp4(variant, 4, 2, 0.75)  

        # here, not enough depth on alternate strand reverse or forward
        variant.INFO['DP4'] = [0,0,4,1]
        assert v.vcf.is_valid_dp4(variant, 4, 2, 0.75) is False
        variant.INFO['DP4'] = [0,0,1,4]
        assert v.vcf.is_valid_dp4(variant, 4, 2, 0.75) is False

        # mimimum ratio must be > 0.75 
        variant.INFO['DP4'] = [25,0,75,75]
        assert v.vcf.is_valid_dp4(variant, 4, 2, 0.75) is True

        variant.INFO['DP4'] = [25,25,75,74]  # just below 0.75 for the alt reverse
        assert v.vcf.is_valid_dp4(variant, 4, 2, 0.75) is False
        variant.INFO['DP4'] = [25,25,74,75] # just below 0.75 for the alt forward
        assert v.vcf.is_valid_dp4(variant, 4, 2, 0.75) is False
    # variant.ALT is equal to "A" 
    validate_variant_alternatate(variant)

    def validate_variant_reference(variant):
        # variant.ALT must be different from "." for this test
        assert str(variant.ALT[0].value).strip() == "."

        # test minimum depth of alternate must be >= 4
        variant.INFO['DP4'] = [2, 2, 0, 0]
        assert v.vcf.is_valid_dp4(variant, 4, 2, 0.75)

        # here, not enough depth on alternate strand reverse or forward
        variant.INFO['DP4'] = [4, 1, 0, 0]
        assert v.vcf.is_valid_dp4(variant, 4, 2, 0.75) is False
        variant.INFO['DP4'] = [1, 4, 0, 0]
        assert v.vcf.is_valid_dp4(variant, 4, 2, 0.75) is False

        # mimimum ratio must be > 0.75
        variant.INFO['DP4'] = [75, 75, 25, 0]
        assert v.vcf.is_valid_dp4(variant, 4, 2, 0.75) is True

        variant.INFO['DP4'] = [75, 74, 25, 25]  # just below 0.75 for the alt reverse
        assert v.vcf.is_valid_dp4(variant, 4, 2, 0.75) is False
        variant.INFO['DP4'] = [74, 75, 25, 25] # just below 0.75 for the alt forward
        assert v.vcf.is_valid_dp4(variant, 4, 2, 0.75) is False

    # variant.ALT is equal to "A"
    variant.ALT[0].value = "."
    validate_variant_reference(variant)


    # Now, let us do the filtering with the vcf_filter method
    v = VCF(data)
    v.vcf.apply_dp4_filter = True
    with TempFile() as fh:
        res = v.vcf.filter_vcf(fh.name)
    assert res ==  {'N': 573, 'filtered': 414, 'unfiltered': 159}

