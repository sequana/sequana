import shutil

import pytest
from easydev import TempFile

from sequana import sequana_data
from sequana.scripts.main import mapping as script

from ... import test_dir


def test_analysis():
    from click.testing import CliRunner

    runner = CliRunner()

    file1 = sequana_data("Hm2_GTGAAA_L005_R1_001.fastq.gz")
    file2 = sequana_data("Hm2_GTGAAA_L005_R2_001.fastq.gz")
    reference = f"{test_dir}/data/fasta/measles.fa"

    from tempfile import TemporaryDirectory

    directory = TemporaryDirectory()
    # shutil.copy(file1, directory.name)
    # shutil.copy(file2, directory.name)
    shutil.copy(reference, directory.name)

    results = runner.invoke(
        script.mapping, ["--file1", file1, "--file2", file2, "--reference", directory.name + "/measles.fa"]
    )
    assert results.exit_code == 0
