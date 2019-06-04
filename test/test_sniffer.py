from sequana.sniffer import sniffer
from sequana import sequana_data


def test_sniffer():
    assert sniffer(sequana_data("test_measles.sam")) == "SAM"
    assert sniffer(sequana_data("test_measles.bam")) == "BAM"
    assert sniffer(sequana_data("test_measles.cram")) == "CRAM"
    assert sniffer(sequana_data("test.fasta")) == "FASTA"
    assert sniffer(sequana_data("test.fastq")) == "FASTQ"
