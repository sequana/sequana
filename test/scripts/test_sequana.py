from easydev import TempFile
from sequana.scripts import main as script
import subprocess

from sequana import sequana_data
filename1 = sequana_data("test.fastq")
filename2 = sequana_data("test_1_1000.fastq.gz")

def test_sequana_app():
    from click.testing import CliRunner
    runner = CliRunner()

    results = runner.invoke(script.fastq, ['--help'])
    # fastq
    results = runner.invoke(script.fastq, [filename1, '--count-reads'])
    assert results.exit_code == 0
    with TempFile() as fout:
        results = runner.invoke(script.fastq, [filename1, '--head', 100, '--output',
fout.name])
        assert results.exit_code == 0


def test_sequana_app_user():

    cmd = "sequana --version"
    status = subprocess.call(cmd.split())
    assert status == 0

    cmd = "sequana fastq --help"
    subprocess.call(cmd.split())

    cmd = ["sequana", "fastq", filename1, "--count-reads"]
    subprocess.call(cmd)
    

    #with TempFile() as fout:
    #    cmd = ["sequana", "fastq", filename1, "--head", '--output', fout.name]
    #    subprocess.call(cmd)


