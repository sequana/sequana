"""Tests for `ribodesigner` package."""

import filecmp
import shutil
from pathlib import Path

import pytest
from click.testing import CliRunner
from sequana.ribodesigner import RiboDesigner
from sequana.scripts.main import ribodesigner

from . import test_dir

resources_dir = Path(test_dir) / "data" / "ribodesigner"


def test_ribodesigner(tmp_path):
    rd = RiboDesigner(
        fasta=resources_dir / "sample.fasta", gff=resources_dir / "sample.gff", output_directory=tmp_path, force=True
    )
    rd.run()
    assert rd.filtered_gff_df.shape == (3, 9)
    assert filecmp.cmp(tmp_path / "probes_sequences.fas", resources_dir / "probes.fas")
    assert filecmp.cmp(tmp_path / "clustered_probes.fas", resources_dir / "probes_clustered.fas")
    assert filecmp.cmp(tmp_path / "clustered_probes.csv", resources_dir / "probes_clustered.csv")


def test_ribodesigner_cli(tmp_path):

    runner = CliRunner()
    result = runner.invoke(
        ribodesigner.ribodesigner,
        [
            str(resources_dir / "sample.fasta"),
            str(resources_dir / "sample.gff"),
            "--output-directory",
            tmp_path / "out_ribodesigner",
        ],
    )
    assert result.exit_code == 0
