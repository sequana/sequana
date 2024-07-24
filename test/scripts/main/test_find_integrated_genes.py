import subprocess

import pytest

from sequana.scripts.main import find_integrated_genes as script

from ... import test_dir


def test_sequana_app(tmpdir):

    from click.testing import CliRunner

    runner = CliRunner()

    # The samplesheet command
    infile = f"{test_dir}/data/bam/measles.fa.sorted.bam"

    print(infile)
    results = runner.invoke(script.find_integrated_genes, ["--help"])
    assert results.exit_code == 0


    results = runner.invoke(script.find_integrated_genes, ["--bam-file", infile, "--name", "ENA|K01711|K01711.1", "--save-reads", "--tag", str(tmpdir) ])

    

    assert results.exit_code == 0

