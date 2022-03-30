from sequana import BED



from . import test_dir

def test_bed():
    bedfile = f"{test_dir}/data/bed/hg38_chr18.bed"
    b = BED(bedfile)
    assert len(list(b.get_transcript_ranges())) == 2777
    assert len(b.get_exons()) == 32640
    assert len(b.get_CDS_exons()) == 24263
    assert len(b) == 2777


def test_wrong_file():
    bedfile = f"{test_dir}/data/bed/test_wrong_bed.bed"
    b = BED(bedfile)
    assert len(b) == 2
    try:
        b.get_exons()
        assert False
    except AssertionError:
        assert True
    except:
        assert False


    # This function uses _get_line that returns {} if a line is incorrectly
    # constructed (12 columns). Previous try/except uses get_exons that simply
    # fails
    b.get_CDS_exons()


    assert b._get_line("1 2 3") == {}




