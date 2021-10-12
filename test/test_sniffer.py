from sequana.sniffer import sniffer
from sequana import sequana_data

from . import test_dir

def test_sniffer():
    assert sniffer(sequana_data("test_measles.sam")) == "SAM"
    assert sniffer(sequana_data("test_measles.bam")) == "BAM"
    assert sniffer(sequana_data("test_measles.cram")) == "CRAM"
    assert sniffer(sequana_data("test.fasta")) == "FASTA"
    assert sniffer(f"{test_dir}/data/fastq/test.fastq") == "FASTQ"
