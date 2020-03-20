from multiqc.utils import config

def load_config():
    my_search_patterns = {
        'sequana_coverage': {'fn': 'sequana_summary_coverage.json'},
        'sequana_pacbio_qc': {'fn': 'sequana_summary*.json'},
        'sequana_quality_control': {'fn': 'summary*.json'},
        'sequana_isoseq_qc': {'fn': 'sequana_summary*.json'},
        'sequana_isoseq': {'fn': 'sequana_summary*.json'},
        'sequana_kraken': {'fn': 'sequana_kraken_summary.json'},
        'sequana_pacbio_amplicon': {'fn': 'sequana_pacbio_amplicon_*.json'},
        'sequana_bamtools_stats': {'fn': 'sequana_bamtools_stats_*.txt'},
    }
    config.update_dict(config.sp, my_search_patterns)

