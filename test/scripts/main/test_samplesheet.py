from easydev import TempFile
from sequana.scripts.main import samplesheet as script
import subprocess
import pytest
from sequana import sequana_data



def test_sequana_app():
    from click.testing import CliRunner
    runner = CliRunner()

    # The samplesheet command
    filename1 = sequana_data("test_expdesign_wrong.csv")
    filename2 = sequana_data("test_expdesign_miseq_illumina_1.csv")
    results = runner.invoke(script.samplesheet, ['--help'])
    results = runner.invoke(script.samplesheet, [filename1, '--check'])
    results = runner.invoke(script.samplesheet, [filename2, '--check'])
    results = runner.invoke(script.samplesheet, [filename2, '--extract-adapters'])


