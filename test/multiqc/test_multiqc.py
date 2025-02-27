from multiqc import report
from multiqc.utils import testing

from . import test_dir


def sequana_data(file):
    return f"{test_dir}/data/{file}"


try:
    from sequana.multiqc import bamtools_stats, coverage, kraken, pacbio_qc

    def test_pacbio():
        # When calling multiqc on the command line, it scans the directory
        # to identify the files to include in the singleton "report";
        # HEre, because we do not use the standalone app, the report.files is empty
        # so we populate it by hand. Moreovoer, the path are altered to look for
        # files in the sequana/resources/testing directory instead of local
        # directory. Because we populate the report.files ourself, we can put
        # whatever name except it the MultiqcModule expects a specific name

        report.reset()
        report.files = {
            "sequana_pacbio_qc": [
                {
                    "filesize": 5913,
                    "fn": sequana_data("summary_pacbio_qc1.json"),
                    "root": ".",
                    "sp_key": "sequana_pacbio_qc",
                },
                {
                    "filesize": 5731,
                    "fn": sequana_data("summary_pacbio_qc2.json"),
                    "root": ".",
                    "sp_key": "sequana_pacbio_qc",
                },
                {
                    "filesize": 5820,
                    "fn": sequana_data("summary_pacbio_qc3.json"),
                    "root": ".",
                    "sp_key": "sequana_pacbio_qc",
                },
            ]
        }
        pacbio_qc.MultiqcModule()

    def test_coverage():
        report.reset()
        report.files = {
            "sequana_coverage": [
                {"fn": sequana_data("summary_coverage1.json"), "root": ".", "sp_key": "sequana_coverage"},
                {"fn": sequana_data("summary_coverage1.json"), "root": ".", "sp_key": "sequana_coverage"},
            ]
        }
        coverage.MultiqcModule()

    def test_sequana_bamtools():
        report.reset()
        report.files = {
            "sequana_bamtools_stats": [
                {"fn": sequana_data("summary_bamtools_stats.txt"), "root": ".", "sp_key": "sequana_bamtools_stats"},
                {"fn": sequana_data("summary_bamtools_stats.txt"), "root": ".", "sp_key": "sequana_bamtools_stats"},
            ]
        }
        bamtools_stats.MultiqcModule()

    def test_kraken():
        report.reset()
        report.files = {
            "sequana_kraken": [
                {"fn": sequana_data("sequana_kraken_summary.json"), "root": ".", "sp_key": "sequana_kraken"},
            ]
        }
        kraken.MultiqcModule()

except (TypeError, AttributeError):
    pass
except:
    raise IOError
