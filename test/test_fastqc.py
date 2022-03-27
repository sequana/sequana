
from . import test_dir

def test_fastqc():

    from sequana.fastqc import FastQC
    f = FastQC()

    filename = f"{test_dir}/data/misc/test_fastqc_report.zip"
    f.read_sample(filename, "test")
    f.plot_sequence_quality()


