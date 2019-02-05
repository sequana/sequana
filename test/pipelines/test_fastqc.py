import os
import json
import tempfile
from sequana import sequana_data
import subprocess

from .common import Pipeline

class FastQCPipeline(Pipeline):
    def __init__(self, wk=None):
        super(FastQCPipeline, self).__init__(wk)
        data = sequana_data("Hm2_GTGAAA_L005_R1_001.fastq.gz")
        input_directory = os.path.dirname(data)
        self.input_pattern = "Hm*gz"
        self.pipeline = "fastqc"

        self.output = "summary.html"
        cmd = [
            "sequana", "--pipeline", self.pipeline, "--input-directory",
input_directory, "--input-pattern", '%s'% self.input_pattern,
            "--working-directory", self.wk, "--force"]

        if "TRAVIS_PYTHON_VERSION" in os.environ:
            cmd += ["--snakemake-jobs", "1"]

        subprocess.check_call(cmd)


    def check(self):
        assert os.path.exists(self.wk + "/index.html")
        assert os.path.exists(self.wk + "/tree.html")
        assert os.path.exists(self.wk + "/multiqc/multiqc_report.html")


def test_fastqc():
    QC = FastQCPipeline()
    try:
        QC.run()
        QC.check()
        QC.clean()
    except:
        QC.clean()
        raise Exception


