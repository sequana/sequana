from sequana.sniffer import sniffer

from . import test_dir


def test_sniffer():
    assert sniffer(f"{test_dir}/data/sam/test_measles.sam") == "SAM"
    assert sniffer(f"{test_dir}/data/bam/test_measles.bam") == "BAM"
    assert sniffer(f"{test_dir}/data/cram/test_measles.cram") == "CRAM"
    assert sniffer(f"{test_dir}/data/fasta/measles.fa") == "FASTA"
    assert sniffer(f"{test_dir}/data/fastq/test.fastq") == "FASTQ"
