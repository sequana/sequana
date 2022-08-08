from easydev import TempFile
from sequana.scripts.main import ribodesigner as script
import subprocess
import pytest


from . import test_dir

def test_sequana_app(tmpdir):
    from click.testing import CliRunner
    runner = CliRunner()

    # gtf fixer
    fasta_filename = f"{test_dir}/../../data/ribodesigner/sample.fas"
    gff_filename = f"{test_dir}/../../data/ribodesigner/sample.gff"

    outdir = tmpdir.mkdir("temp")


    with TempFile() as fout:
        results = runner.invoke(script.ribodesigner,
            [fasta_filename, gff_filename, "--output-directory", outdir])

