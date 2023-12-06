import pytest

from sequana.scripts.main import fasta as script


def test_main():
    from click.testing import CliRunner

    runner = CliRunner()

    results = runner.invoke(script.fasta, ["--help"])
    assert results.exit_code == 0
