from easydev import TempFile
from sequana.scripts.main import gtf_fixer as script
import subprocess
import pytest
from sequana import sequana_data



def test_sequana_app():
    from click.testing import CliRunner
    runner = CliRunner()

    # gtf fixer
    filename = sequana_data("test_gtf_fixer.gtf")
    with TempFile() as fout:
        results = runner.invoke(script.gtf_fixer, 
            ['--input', filename, "--output", fout.name])

