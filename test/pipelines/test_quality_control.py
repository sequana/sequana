import os
import json
import tempfile
from sequana import sequana_data
import subprocess

from .common import Pipeline

class QualityPipeline(Pipeline):
    def __init__(self, wk=None):
        super(QualityPipeline, self).__init__(wk)
        data = sequana_data("Hm2_GTGAAA_L005_R1_001.fastq.gz")
        input_directory = os.path.dirname(data)
        self.input_pattern = "Hm*gz"
        self.pipeline = "quality_control"

        self.output = self.wk + "/Hm2_GTGAAA_L005/report_qc_Hm2_GTGAAA_L005/summary.json"
        cmd = [
            "sequana", "--pipeline", self.pipeline,
            "--input-pattern", '%s'% self.input_pattern,
            "--input-directory", '%s'% input_directory,
            "--working-directory", self.wk,
            "--adapters", "Nextera", "--force"]

        if "TRAVIS_PYTHON_VERSION" in os.environ:
            cmd += ["--snakemake-jobs", "1"]

        subprocess.check_call(cmd)


    def check(self):
        if os.path.exists(self.output):
            data = json.load(open(self.output))
        assert data['project'] == "Hm2_GTGAAA_L005" 
        #assert data == truth

        assert data["cutadapt_json"] == {"Number of reads": {'Pairs kept': '1,316',
          'Pairs too short': '175',
          'Read1 with adapters': '104',
          'Read2 with adapters': '142',
          'Total paired reads': '1,491'},
         'percent': {'Pairs kept': '(88.3%)',
          'Pairs too short': '(11.7%)',
          'Read1 with adapters': '(7.0%)',
          'Read2 with adapters': '(9.5%)',
          'Total paired reads': '(100%)'}}

        assert data['phix_section'] == {'R1_mapped': 8,
             'R1_unmapped': 1491,
             'R2_mapped': 8,
             'R2_unmapped': 1491,
             'contamination': 0.5336891260840559,
             'duplicated': 0,
             'mode': 'pe',
             'unpaired': 0}

        return data


def test_quality_control():
    QC = QualityPipeline()
    try:
        QC.run()
        QC.check()
        QC.clean()
    except:
        QC.clean()
        raise Exception


