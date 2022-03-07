from easydev import TempFile
from sequana.scripts.main import lane_merging as script
import subprocess
import pytest
from sequana import sequana_data


from . import test_dir

def test_sequana_app(tmpdir):
    from click.testing import CliRunner
    runner = CliRunner()

    # fastq command
    file1 = f"{test_dir}/data/test_L001_R1_001.fastq.gz"
    file2 = f"{test_dir}/data/test_L002_R1_001.fastq.gz"


    directory_data = tmpdir.mkdir("temp")

    results = runner.invoke(script.lane_merging, ['--help'])
    results = runner.invoke(script.lane_merging, ["--lanes", "1", "2", 
            "--pattern", f"{test_dir}/test*L00*gz", "--output-directory", directory_data.__str__() ])
    assert results.exit_code == 0


