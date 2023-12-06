import pytest

from sequana.scripts.main import rnadiff as script


def test_main():
    from click.testing import CliRunner

    runner = CliRunner()

    results = runner.invoke(script.rnadiff, ["--help"])
    assert results.exit_code == 0
