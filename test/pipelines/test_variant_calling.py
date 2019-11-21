from sequana import SequanaConfig, sequana_data
from easydev import shellcmd
import subprocess
import json
import os
import tempfile
from .common import Pipeline


class VariantCallingPipeline(Pipeline):

    def __init__(self, wk=None):
        super(VariantCallingPipeline, self).__init__(wk=wk)
        # Define the data
        data = sequana_data("Hm2_GTGAAA_L005_R1_001.fastq.gz")
        input_directory = os.path.dirname(data)
        self.input_pattern = "Hm*gz"
        self.pipeline = "variant_calling"

        # Define the project and config file
        subprocess.check_call([
            "sequana", "--pipeline", self.pipeline,
            "--input-pattern", '%s' % self.input_pattern,
            "--input-directory", '%s' % input_directory,
            "--working-directory", self.wk, "--force"
            ])


        cmd = ["sequana", "--pipeline", self.pipeline,
             "--input-pattern", '%s'% self.input_pattern,
            "--input-directory", '%s' % input_directory,
             "--working-directory", self.wk, "--force"]

        if "TRAVIS_PYTHON_VERSION" in os.environ:
             cmd += ["--snakemake-jobs", "1"]

        subprocess.check_call(cmd)


        # Add reference in the config
        cfg = SequanaConfig(self.wk + "/config.yaml")
        # We added a TTTT in position 5881
        cfg._yaml_code['bwa_mem_ref']['reference'] = sequana_data("measles.fa")
        cfg.save(self.wk + '/config.yaml')

    def check(self):
        from sequana.freebayes_vcf_filter import VCF_freebayes, Variant
        vcf = VCF_freebayes(self.wk +
            "/report_vc_Hm2_GTGAAA_L005/outputs/Hm2_GTGAAA_L005.raw.vcf")
        vcf.rewind()
        vv = [Variant(v)._resume for v in vcf]
        event276 = {'alternative': 'G',
                'chr': 'ENA|K01711|K01711.1',
                'depth': 5,
                'freebayes_score': 126.901,
                'frequency': '1.00',
                'position': '276',
                'reference': 'C',
                'strand_balance': '0.40'}
        event = [v for v in vv if v['position'] =="276"][0]
        assert len(vv) in (85,93) # 85 in freebayes 1.0 and 93 in freebayes 1.2
        for k in ["depth", "chr", "frequency", "position", "reference",
                  "alternative", "strand_balance"]:
            assert event[k] == event276[k]

def test_variant_calling():
    QC = VariantCallingPipeline()
    try:
        QC.run()
        QC.check()
        QC.clean()
    except Exception as err:
        QC.clean()
        raise(err)
