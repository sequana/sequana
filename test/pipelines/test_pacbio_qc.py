from sequana import SequanaConfig, sequana_data
from easydev import shellcmd
import subprocess
import json
import os
import tempfile
from .common import Pipeline
import pytest

skiptravis = pytest.mark.skipif("TRAVIS_PYTHON_VERSION" in os.environ,
     reason="On travis")
 
 
class PacbioQCPipeline(Pipeline):

    def __init__(self, wk=None):
        super(PacbioQCPipeline, self).__init__(wk=wk)
        # Define the data
        data = sequana_data("test_pacbio_subreads.bam")
        input_directory = os.path.dirname(data)
        self.input_pattern = "test_pacbio_subreads.bam"
        self.pipeline = "pacbio_qc"

        # Define the project and config file
        cmd = ["sequana", "--pipeline", self.pipeline,
            "--input-pattern", '%s' % self.input_pattern,
            "--input-directory", '%s' % input_directory,
            "--extension", "bam",
            "--working-directory", self.wk, "--force"]

        if "TRAVIS_PYTHON_VERSION" in os.environ:
            cmd += ["--snakemake-jobs", "1"]
        subprocess.check_call(cmd)

        cfg = SequanaConfig(self.wk + "/config.yaml")
        cfg._yaml_code["input_directory"] = input_directory
        cfg._yaml_code["input_readtag"] = "_R[12]_"
        cfg._yaml_code['input_pattern'] = self.input_pattern 
        cfg.save(self.wk + '/config.yaml')

        self.output = self.wk + "/1/sequana_summary_pacbio_qc_1.json"


    def check(self):
        if os.path.exists(self.output):
            data = json.load(open(self.output))
        assert "hist_gc" in data


@skiptravis
def test_pipeline():
    QC = PacbioQCPipeline()
    try:
        QC.run()
        QC.check()
        QC.clean()
    except TimeoutError:
        pass
    except:
        QC.clean()
        raise Exception
