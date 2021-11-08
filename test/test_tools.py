from sequana.tools import bam_to_mapped_unmapped_fastq, reverse_complement, StatsBAM2Mapped
from sequana.tools import fast_gc_content
from sequana import sequana_data
from sequana.tools import bam_get_paired_distance, GZLineCounter, PairedFastQ
from sequana.tools import genbank_features_parser

def test_StatsBAM2Mapped():
    data = sequana_data("test.bam", "testing")
    res = StatsBAM2Mapped(data)
    res.to_html()

def test_fast_gc_content():
    assert fast_gc_content("GCGCGCGCGC") == 1
    assert fast_gc_content("AAAAAAAAAA") == 0
    assert fast_gc_content("gcgcAAAA") == 0.5
    assert fast_gc_content("gcgcxXXXAAAA") == 0.5

def test_bam2fastq():
    data = sequana_data("test.bam", "testing")
    res = bam_to_mapped_unmapped_fastq(data)


def test_reverse_complement():
    assert reverse_complement("AACCGGTTA") == 'TAACCGGTT'


def test_reverse():
    from sequana.tools import reverse
    assert reverse("AACCGG") == 'GGCCAA'


def test_distance():
    data = sequana_data("test.bam", "testing")
    distances = bam_get_paired_distance(data)


def test_gc_content():
    from sequana.tools import gc_content
    data = sequana_data('measles.fa', "testing")
    gc_content(data, 10)['chr1']
    gc_content(data, 101, circular=True)['chr1']


def test_genbank_features_parser():
    data = sequana_data("JB409847.gbk")
    genbank_features_parser(data)


def test_gzlinecounter():
    assert len(GZLineCounter(sequana_data("test.fastq.gz"))) == 1000


def test_paired_file():
    f1 = sequana_data("test.fastq.gz")
    f2 = sequana_data("test.fastq.gz")
    assert PairedFastQ(f1,f2).is_synchronised()
