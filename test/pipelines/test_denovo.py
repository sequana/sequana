import os
import tempfile
from sequana import sequana_data
import subprocess
from sequana import SequanaConfig
from .common import Pipeline


class DenovoPipeline(Pipeline):
    def __init__(self, wk=None):
        super(DenovoPipeline, self).__init__(wk)
        data = sequana_data("Hm2_GTGAAA_L005_R1_001.fastq.gz")
        input_directory = os.path.dirname(data)
        self.input_pattern = input_directory + "/Hm*gz"
        self.pipeline = "denovo_assembly"

        #self.output = self.wk + "/Hm2_GTGAAA_L005/report_qc_Hm2_GTGAAA_L005/summary.json"
        subprocess.check_call([
            "sequana", "--pipeline", self.pipeline,
            "--input-pattern", '%s'% self.input_pattern,
            "--working-directory", self.wk,
            "--force"])

        # Add reference in the config
        cfg = SequanaConfig(self.wk + "/config.yaml")
        # We added a TTTT in position 5881
        cfg._yaml_code['digital_normalisation']['max_memory_usage'] = 1e9
        cfg.save(self.wk + '/config.yaml')


    def check(self):
        assert os.path.exists(self.wk + "/quast")


def test_denovo():
    QC = DenovoPipeline()
    try:
        QC.run()
        QC.check()
        QC.clean()
    except:
        QC.clean()
        assert Exception



