from sequana import sequana_data


def test_fastqc():

    from sequana.fastqc import FastQC
    f = FastQC()

    filename = sequana_data("test_fastqc_report.zip")
    f.read_sample(filename, "test")
    f.plot_sequence_quality()


