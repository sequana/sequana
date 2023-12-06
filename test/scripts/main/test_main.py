import pytest

from sequana.scripts.main import main as script


def test_main():
    from click.testing import CliRunner

    runner = CliRunner()

    results = runner.invoke(script.main, ["--help"])
    assert results.exit_code == 0
