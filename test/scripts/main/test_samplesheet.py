import subprocess

import pytest
from easydev import TempFile

from sequana.scripts.main import samplesheet as script

from ... import test_dir


def test_sequana_app(tmpdir):
    from click.testing import CliRunner

    runner = CliRunner()

    # The samplesheet command
    filename1 = f"{test_dir}/data/csv/test_expdesign_wrong.csv"
    filename2 = f"{test_dir}/data/csv/test_expdesign_miseq_illumina_1.csv"
    results = runner.invoke(script.samplesheet, ["--help"])
    assert results.exit_code == 0

    results = runner.invoke(script.samplesheet, [filename1, "--check"])
    assert results.exit_code == 0

    results = runner.invoke(script.samplesheet, [filename2, "--check"])
    assert results.exit_code == 0

    output = tmpdir.mkdir("test").join("test.csv")
    results = runner.invoke(script.samplesheet, [filename2, "--quick-fix", "--output", output])
    assert results.exit_code == 0

    results = runner.invoke(script.samplesheet, [filename2, "--extract-adapters"])
    assert results.exit_code == 0
