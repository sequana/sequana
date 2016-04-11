from sequana.bamtools import BAM, Alignment, SAMFlags
from sequana.report_bam import BAMReport
from sequana import sequana_data
from easydev import TempFile


datatest = sequana_data("test.bam", "testing")


def test_bam():
    s = BAM(datatest)
    assert len(s) == 1000

    assert len(list(s.iter_unmapped_reads())) == 2
    s.reset()
    assert len(list(s.iter_mapped_reads())) == 998
    s.reset()

    assert s.get_read_names()

    s.get_stats()
    s.get_full_stats_as_df()

    with TempFile(suffix='.png') as fh:
        s.plot_bar_flags(filename=fh.name, logy=True)
        s.plot_bar_flags(filename=fh.name)

    with TempFile(suffix='.png') as fh:
        s.plot_bar_mapq(filename=fh.name)

    with TempFile() as fh:
        s.to_fastq(fh.name)
        from sequana import FastQ
        ff = FastQ(fh.name)
        len(ff) == len(s)
    


def test_alignment():
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


def test_bamreport():
    b = BAM(datatest)
    r = BAMReport()
    r.set_data(b)
    r.create_report()
