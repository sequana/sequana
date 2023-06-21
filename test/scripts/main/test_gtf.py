from sequana.scripts.main import gtf_fixer as script
import subprocess
import pytest


from . import test_dir

def test_sequana_app(tmpdir):

    fout = tmpdir.join("out.gtf")
    from click.testing import CliRunner
    runner = CliRunner()

    # gtf fixer
    filename = f"{test_dir}/data/test_gtf_fixer.gtf"
    results = runner.invoke(script.gtf_fixer, ['--input', filename, "--output", fout])

