from sequana.tools import (
    GZLineCounter,
    PairedFastQ,
    StatsBAM2Mapped,
    bam_get_paired_distance,
    bam_to_mapped_unmapped_fastq,
    entropy,
    fast_gc_content,
    gc_content,
    genbank_features_parser,
    reverse,
    reverse_complement,
)

from . import test_dir


def test_StatsBAM2Mapped(tmpdir):

    data = f"{test_dir}/data/bam/test.bam"
    res = StatsBAM2Mapped(data)
    res.to_html()

    json_filename = tmpdir.join("test.json")
    res.to_json(json_filename)


def test_entropy():
    entropy("ACGT")


def test_fast_gc_content():
    assert fast_gc_content("GCGCGCGCGC") == 1
    assert fast_gc_content("AAAAAAAAAA") == 0
    assert fast_gc_content("gcgcAAAA") == 0.5
    assert fast_gc_content("gcgcxXXXAAAA") == 0.5


def test_bam2fastq():
    data = f"{test_dir}/data/bam/test.bam"
    res = bam_to_mapped_unmapped_fastq(data)


def test_reverse_complement():
    assert reverse_complement("AACCGGTTA") == "TAACCGGTT"


def test_reverse():
    assert reverse("AACCGG") == "GGCCAA"


def test_distance():
    data = f"{test_dir}/data/bam/test.bam"
    distances = bam_get_paired_distance(data)


def test_gc_content():
    data = f"{test_dir}/data/fasta/measles.fa"
    gc_content(data, 10)["chr1"]
    gc_content(data, 101, circular=True)["chr1"]


def test_genbank_features_parser():
    data = f"{test_dir}/data/genbank/JB409847.gbk"
    genbank_features_parser(data)


def test_gzlinecounter():
    filename = f"{test_dir}/data/fastq/test.fastq.gz"
    g = GZLineCounter(filename)
    assert len(g) == 1000

    g.use_zcat = False
    assert len(g) == 1000


def test_paired_file():
    f1 = f"{test_dir}/data/fastq/test.fastq.gz"
    f2 = f"{test_dir}/data/fastq/test.fastq.gz"
    assert PairedFastQ(f1, f2).is_synchronised()

    f1 = f"{test_dir}/data/fastq/test_R1.fastq"
    f2 = f"{test_dir}/data/fastq/test_R2.fastq"
    assert not PairedFastQ(f1, f2).is_synchronised()
