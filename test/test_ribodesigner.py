#!/usr/bin/env python

"""Tests for `ribodesigner` package."""

import pytest

from click.testing import CliRunner

from sequana import ribodesigner as rd
from sequana.scripts import main

from . import test_dir
from pathlib import Path
import shutil

resources_dir = Path(test_dir) / "data" / "ribodesigner"


@pytest.fixture(scope="function")
def setup_ribo(tmpdir):
    shutil.copytree(resources_dir, tmpdir, dirs_exist_ok=True)


def test_get_rna_pos_from_gff(setup_ribo, tmpdir):
    df = rd.get_rna_pos_from_gff(tmpdir / "sample.gff")
    print(df)
    assert df.shape == (3, 9)


def test_extract_regions_from_fasta(setup_ribo, tmpdir):
    df = rd.get_rna_pos_from_gff(tmpdir / "sample.gff")
    rd.extract_regions_from_fasta(tmpdir / "sample.fasta", df, tmpdir / "test_regions.fas")
    with open(tmpdir / "test_regions.fas") as f1:
        with open(tmpdir / "regions.fas") as f2:
            assert list(f1) == list(f2)


def test_reverse_complement():

    sequence = "CGTAGTGCTAGTGCTAGTGTCGATGTCGATGTCGATGTCGATCTAGTC"
    rev_comp_sequence = "GACTAGATCGACATCGACATCGACATCGACACTAGCACTAGCACTACG"

    assert rd.reverse_complement(sequence) == rev_comp_sequence


def test_get_probes(setup_ribo, tmpdir):

    rd.get_probes(tmpdir / "regions.fas", tmpdir / "test_probes.fas")
    with open(tmpdir / "test_probes.fas") as f1:
        with open(tmpdir / "probes.fas") as f2:
            assert list(f1) == list(f2)


def test_cluster_probes(setup_ribo, tmpdir):

    rd.cluster_probes(tmpdir / "probes.fas", tmpdir / "test_probes_clustered.fas")
    with open(tmpdir / "test_probes_clustered.fas") as f1:
        with open(tmpdir / "probes_clustered.fas") as f2:
            assert list(f1) == list(f2)


def test_fasta_to_csv(setup_ribo, tmpdir):

    rd.fasta_to_csv(tmpdir / "probes_clustered.fas", tmpdir / "test_probes_clustered.csv")
    with open(tmpdir / "test_probes_clustered.csv") as f1:
        with open(tmpdir / "probes_clustered.csv") as f2:
            assert list(f1) == list(f2)


def test_ribodesigner_cli(setup_ribo, tmpdir):

    runner = CliRunner()
    result = runner.invoke(
        main.ribodesigner,
        [str(tmpdir / "sample.fasta"), str(tmpdir / "sample.gff"), "--outdir", tmpdir / "out_ribodesigner"],
    )
    print(result.output)
    assert result.exit_code == 0
