from easydev import TempFile
from sequana.scripts.main import fastq as script
import subprocess
import pytest
from sequana import sequana_data


from . import test_dir

def test_sequana_app():
    from click.testing import CliRunner
    runner = CliRunner()

    # fastq command
    filename1 = f"{test_dir}/../../data/fastq/test.fastq.gz"
    results = runner.invoke(script.fastq, ['--help'])
    results = runner.invoke(script.fastq, [filename1, '--count-reads'])
    assert results.exit_code == 0
    with TempFile() as fout:
        results = runner.invoke(script.fastq, [filename1, 
            '--head', 100, '--output', fout.name])
        assert results.exit_code == 0


