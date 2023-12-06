import pytest

from sequana.scripts.main import biomart as script


def test_main():
    from click.testing import CliRunner

    runner = CliRunner()

    results = runner.invoke(script.biomart, ["--help"])
    assert results.exit_code == 0
