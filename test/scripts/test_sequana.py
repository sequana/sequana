from easydev import TempFile
from sequana.scripts import main as script
import subprocess

from sequana import sequana_data

def test_sequana_app():
    from click.testing import CliRunner
    runner = CliRunner()

    # fastq command
    filename1 = sequana_data("test.fastq")
    filename2 = sequana_data("test_1_1000.fastq.gz")
    results = runner.invoke(script.fastq, ['--help'])
    results = runner.invoke(script.fastq, [filename1, '--count-reads'])
    assert results.exit_code == 0
    with TempFile() as fout:
        results = runner.invoke(script.fastq, [filename1, 
            '--head', 100, '--output', fout.name])
        assert results.exit_code == 0

    # The samplesheet command

    filename1 = sequana_data("test_expdesign_wrong.csv")
    filename2 = sequana_data("test_expdesign_miseq_illumina_1.csv")
    results = runner.invoke(script.samplesheet, ['--help'])
    results = runner.invoke(script.samplesheet, [filename1, '--check'])
    results = runner.invoke(script.samplesheet, [filename2, '--check'])
    results = runner.invoke(script.samplesheet, [filename2, '--extract-adapters'])

    # gtf fixer
    filename = sequana_data("test_gtf_fixer.gtf")
    with TempFile() as fout:
        results = runner.invoke(script.gtf_fixer, 
            ['--input', filename, "--output", fout.name])
    


def test_sequana_app_user():
    filename1 = sequana_data("test.fastq")
    filename2 = sequana_data("test_1_1000.fastq.gz")

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


