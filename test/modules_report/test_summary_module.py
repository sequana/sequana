import json
from sequana import sequana_data
from sequana.modules_report.summary import SummaryModule
from sequana.utils import config
from . import test_dir

sharedir = f"{test_dir}/data"

def test_summary_module(tmpdir):
    directory = tmpdir.mkdir('test_variant_calling_module')
    config.output_dir = str(directory)
    config.sample_name = 'JB409847'
    summary_dict = {'tool': 'sequana_summary',
                    'inputs': [
                        sequana_data('Hm2_GTGAAA_L005_R1_001.fastq.gz'),
                        sequana_data('Hm2_GTGAAA_L005_R2_001.fastq.gz')],
                    'outputs': [f"{sharedir}/JB409847.vcf"],
                    'html': [f"{sharedir}/JB409847.vcf"],
                    'rulegraph': f"{sharedir}/test_summary_module.svg",
                    'requirements': f"{sharedir}/test_summary_module.svg",
                    'snakefile': f"{sharedir}/test_summary_module.svg",
                    'config': f"{sharedir}/test_summary_module.svg",
                    'name': 'JB409847'}
    SummaryModule(summary_dict)
