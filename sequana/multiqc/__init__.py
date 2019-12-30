from __future__ import absolute_import
from multiqc import config

# Add search patterns and config options for the things that are used in
# MultiQC_sequana
def multiqc_sequana_config():
    """ Set up MultiQC config defaults for this package """
    sequana_search_patterns = {
        'sequana/pacbio_qc': {
            'fn': 'summary_*.json'
        },
        'sequana/quality_control': {
            'fn': 'summary.json'
        },
        'sequana/coverage': {
            'fn': 'sequana_summary_coverage.json',
        },
        'sequana/isoseq': {
            'fn': 'sequana_summary_isoseq.json',
        },
        'sequana/isoseq_qc': {
            'fn': 'sequana_summary_isoseq_qc.json',
        },
        'sequana/pacbio_amplicon': {
            'fn': 'sequana_pacbio_amplicon*.json',
        },
    }
    config.update_dict(config.sp, sequana_search_patterns)
    #config.fn_clean_exts.append({'type': 'regex', 'pattern': 'summary_*.*'})
    config.update_dict(config.table_columns_visible, {
        'FastQC': {
            'percent_duplicates': False,
            'total_sequences': False,
        },
        'Samtools Stats': {
            'non-primary_alignments': False,
            'reads_mapped': False,
            'reads_mapped_percent': False,
            'raw_total_sequences': False,
        }
    })

    # config.log_filesize_limit =  500000
    # in sequana_coverage, large html files are created.
    # they can be ignored
    config.fn_ignore_files = ['*html']


    #config.update_dict(config.log_filesize_limit, 500000)
