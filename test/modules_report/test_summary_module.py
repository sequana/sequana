import json

from sequana import sequana_data
from sequana.modules_report.summary import SummaryModule
from sequana.utils import config


def test_summary_module(tmpdir):
    directory = tmpdir.mkdir('test_variant_calling_module')
    config.output_dir = str(directory)
    config.sample_name = 'JB409847'
    summary_dict = {'tool': 'sequana_summary',
                    'inputs': [
                        sequana_data('Hm2_GTGAAA_L005_R1_001.fastq.gz'),
                        sequana_data('Hm2_GTGAAA_L005_R2_001.fastq.gz')],
                    'outputs': [sequana_data('JB409847.vcf')],
                    'html': [sequana_data('JB409847.vcf')],
                    'rulegraph': sequana_data('test_summary_module.svg'),
                    'requirements': sequana_data('test_summary_module.svg'),
                    'snakefile': sequana_data('test_summary_module.svg'),
                    'config': sequana_data('test_summary_module.svg'),
                    'name': 'JB409847'}
    SummaryModule(summary_dict)
