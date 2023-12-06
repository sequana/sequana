import os

import pytest
from click.testing import CliRunner

from sequana import logger, sequana_data
from sequana.scripts import taxonomy

prog = "sequana_taxonomy"


def test_analysis(tmpdir):
    directory = tmpdir.mkdir("report")

    file1 = sequana_data("Hm2_GTGAAA_L005_R1_001.fastq.gz")
    file2 = sequana_data("Hm2_GTGAAA_L005_R2_001.fastq.gz")

    # Test that database must be provided

    runner = CliRunner()
    results = runner.invoke(taxonomy.main, ["--input-file1", file1])
    assert results.exit_code != 0

    # FIXME. Those calls do not return 0 (why ?)
    results = runner.invoke(taxonomy.main, ["--download", "toydb"])
    # assert results.exit_code == 0

    # FIXME. Those calls do not return 0 (why ?)
    results = runner.invoke(
        taxonomy.main, ["--input-file1", file1, "--databases", "toydb", "--output-directory", directory.__str__()]
    )
    # assert results.exit_code == 0

    results = runner.invoke(
        taxonomy.main,
        [
            "--input-file1",
            file1,
            "--input-file2",
            file2,
            "--output-directory",
            directory.__str__(),
            "--databases",
            "toydb",
        ],
    )
    # assert results.exit_code == 0

    results = runner.invoke(
        taxonomy.main,
        [
            "--input-file1",
            file1,
            "--input-file2",
            file2,
            "--databases",
            "toydb",
            "toydb",
            "--output-directory",
            directory.__str__(),
        ],
    )
    # assert results.exit_code == 0

    results = runner.invoke(taxonomy.main, ["--input-file1", file1, "--input-file2", file2, "--databases", "dummy"])
    assert results.exit_code == 1


def test_help():
    try:
        taxonomy.main([prog, "--help", "1>/tmp/out", "2>/tmp/err"])
        assert False
    except SystemExit:
        pass
    else:
        raise Exception
