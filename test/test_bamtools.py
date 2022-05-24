from sequana.bamtools import SAM, BAM, CRAM,  Alignment, SAMFlags
from sequana.bamtools import is_bam, is_cram, is_sam
from sequana.modules_report.bamqc import BAMQCModule
from easydev import TempFile
import pytest
import os

from . import test_dir

def test_is_sam_bam():

    datatest = f"{test_dir}/data/sam/test_measles.sam"
    assert is_sam(datatest) is True

    datatest = f"{test_dir}/data/bam/test_measles.bam"
    assert is_bam(datatest) is True

def test_is_cram():
    datatest = f"{test_dir}/data/cram/test_measles.cram"
    assert is_cram(datatest) is True


def test_sam(tmpdir):
    datatest = f"{test_dir}/data/sam/test.sam"
    s = SAM(datatest)
    assert len(s) == 432
    assert s.is_sorted is True
    assert s.is_paired is True
    df = s.get_df_concordance(max_align=100)
    s.hist_soft_clipping()


def test_bam(tmpdir):
    datatest = f"{test_dir}/data/bam/test.bam"
    s = BAM(datatest)
    assert len(s) == 1000
    assert s.is_sorted is True
    df = s.get_df_concordance()
    assert s.is_paired is True
    assert int(df.length.sum()) == 67938
    assert int(df.M.sum()) == 67788

    df = s.get_df()

    # call this here before other computations on purpose
    with TempFile(suffix=".json") as fh:
        s.bam_analysis_to_json(fh.name)

    assert s.get_read_names()
    s.get_mapped_read_length()

    s.get_stats()
    s.get_stats_full()
    s.get_samtools_stats_as_df()

    with TempFile() as fh:
        s.to_fastq(fh.name)
        from sequana import FastQ
        ff = FastQ(fh.name)
        len(ff) == len(s)

    # plotting
    with TempFile(suffix='.png') as fh:
        s.plot_bar_flags(filename=fh.name, logy=True)
        s.plot_bar_flags(filename=fh.name)

    with TempFile(suffix='.png') as fh:
        s.plot_bar_mapq(filename=fh.name)

    s.get_gc_content()
    s.get_length_count()
    s.plot_gc_content()
    s.boxplot_qualities()
    s.boxplot_qualities(max_sample=50)
    try:
        s.plot_gc_content(bins=[1,2,10])
        assert False
    except:
        assert True

@pytest.mark.xfail(reason="cram flaky")
def test_cram():
    datatest = f"{test_dir}/data/cram/test_measles.cram"
    s = CRAM(datatest)
    assert s.summary == {'flags': {77: 6, 83: 14, 99: 10, 141: 6, 147: 10, 163: 14},
         'mapq': {0: 12, 60: 48},
         'mean_quality': 33.666171617161723,  
         'read_length': {79: 2, 81: 1, 93: 1, 101: 44}}

def test_bam_others():
    datatest = f"{test_dir}/data/bam/measles.fa.sorted.bam"
    b = BAM(datatest)
    assert len(b) == 2998

    # plot_reaqd_length and data
    X, Y = b._get_read_length()
    assert sum(Y) == 2623
    b.plot_read_length()
    b.hist_coverage()
    b.plot_coverage()
    b.boxplot_qualities()
    b.plot_indel_dist()


def test_alignment():
    datatest = f"{test_dir}/data/bam/test.bam"
    s = BAM(datatest)
    # no need to call reset but does not harm and reminds us that it shoudl be
    # used in general to make sure we start at the beginning of the iterator.
    s.reset()
    a = Alignment(next(s))
    a.as_dict()


def test_samflags():
    sf = SAMFlags(4095)
    print(sf)
    sf.get_meaning()
    sf.get_flags()


def test_bamreport(tmpdir):
    datatest = f"{test_dir}/data/bam/test.bam"
    directory = tmpdir.mkdir("bam")
    from sequana.utils import config
    config.output_dir = directory.__str__()
    r = BAMQCModule(datatest, "bam.html")


def test_cs():
    from sequana.bamtools import CS
    assert  CS('-a:6-g:14+g:2+c:9*ac:10-a:13-a') ==  {'D': 4, 'I': 2, 'M': 54, 'S': 1}
    assert  CS('-a:6-g:14+g:2+c:9*ac:10-a:13-a:2*ac') ==  {'D': 4, 'I': 2, 'M': 56, 'S': 2}
    assert  CS('-a:6-g:14+g:2+c:9*ac:10-a:13-a:1+g') ==  {'D': 4, 'I': 3, 'M': 55, 'S': 1}


def test_cs_in_bam():
    b = BAM(f"{test_dir}/data/bam/test_CS_tiny.bam")
    assert  b.summary == {
        'flags': {0: 2, 16: 2},
         'mapq': {60: 4},
         'mean_quality': 0.0,
         'read_length': {1772: 1, 10779: 1, 13726: 1, 20480: 1}}
    df = b.get_df_concordance()
    import math
    del df['rname']
    assert math.floor(df.sum().sum()) == 103813  # exact is 103769.5600734975


def test_insert_size():

    d1 = f"{test_dir}/data/sam/test_measles.sam"
    d3 = f"{test_dir}/data/bam/test.bam"
    d4 = f"{test_dir}/data/bam/test_CS_tiny.bam"


    b1  = BAM(d1)
    # test max_entries
    assert len(b1._get_insert_size_data(10)) == 7
    b1.get_estimate_insert_size(100)
    b1.get_estimate_insert_size()
    b1.plot_insert_size()

    #b2  = BAM(d2)
    #b2.get_estimate_insert_size()
    
    b3  = BAM(d3)
    b3.get_estimate_insert_size()
    
    b4  = BAM(d4)
    assert b4.get_estimate_insert_size() == 0


def test_strandness():
    b = BAM(f"{test_dir}/data/bam/test_hg38_chr18.bam")
    res = b.infer_strandness(f"{test_dir}/data/bed/hg38_chr18.bed", 200000)
    assert res[0] == 'Paired-end'
    assert res[1]> 0.94
    assert res[2]<0.06
    assert res[3] < 0.0011



def test_mRNA_inner_distance():
    b = BAM(f"{test_dir}/data/bam/test_hg38_chr18.bam")
    df = b.mRNA_inner_distance(f"{test_dir}/data/bed/hg38_chr18.bed")
    # Total read pairs  used 382
    # mean insert size: 88.3975155279503
    assert df[0]['val'].mean() >1436 and df[0]['val'].mean()<1437

