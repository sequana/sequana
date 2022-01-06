"""Tests for `ribodesigner` package."""

import filecmp
import shutil
from pathlib import Path

import pytest
from click.testing import CliRunner
from sequana import ribodesigner as rd
from sequana.scripts import main

from . import test_dir

resources_dir = Path(test_dir) / "data" / "ribodesigner"


def test_get_rna_pos_from_gff(tmp_path):
    df = rd.get_rna_pos_from_gff(
        resources_dir / "sample.gff",
        tmp_path / "sample_filtered.gff",
        resources_dir / "sample.fasta",
        tmp_path / "test_regions.fas",
    )
    assert df.shape == (3, 9)


def test_get_probes(tmp_path):

    rd.get_probes(resources_dir / "regions.fas", tmp_path / "test_probes.fas")
    assert filecmp.cmp(tmp_path / "test_probes.fas", resources_dir / "probes.fas")


def test_cluster_probes(tmp_path):

    rd.cluster_probes(resources_dir / "probes.fas", tmp_path / "test_probes_clustered.fas")
    assert filecmp.cmp(tmp_path / "test_probes_clustered.fas", resources_dir / "probes_clustered.fas")


def test_fasta_to_csv(tmp_path):

    rd.fasta_to_csv(resources_dir / "probes_clustered.fas", tmp_path / "test_probes_clustered.csv")
    assert filecmp.cmp(tmp_path / "test_probes_clustered.csv", resources_dir / "probes_clustered.csv")


def test_ribodesigner_cli(tmp_path):

    runner = CliRunner()
    result = runner.invoke(
        main.ribodesigner,
        [
            str(resources_dir / "sample.fasta"),
            str(resources_dir / "sample.gff"),
            "--output-directory",
            tmp_path / "out_ribodesigner",
        ],
    )
    assert result.exit_code == 0
