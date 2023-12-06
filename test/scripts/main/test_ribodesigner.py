import pytest

from sequana.scripts.main import ribodesigner as script


def test_main():
    from click.testing import CliRunner

    runner = CliRunner()

    results = runner.invoke(script.ribodesigner, ["--help"])
    assert results.exit_code == 0
