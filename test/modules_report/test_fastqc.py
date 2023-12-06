from sequana.modules_report import fastqc


def test_fastqc(tmpdir):
    directory = tmpdir.mkdir("fastqc")
    from sequana.utils import config

    config.output_dir = directory.__str__()
    ff = fastqc.FastQCModule()
