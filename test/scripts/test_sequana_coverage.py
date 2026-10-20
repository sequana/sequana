import os

import pandas as pd
import pytest
from click.testing import CliRunner

from sequana.scripts import coverage

prog = "sequana_coverage"

from .. import test_dir

bedfile = f"{test_dir}/data/bed/JB409847.bed"
genbank = f"{test_dir}/data/genbank/JB409847.gbk"
fastafile = f"{test_dir}/data/fasta/JB409847.fasta"


def test_run():

    runner = CliRunner()

    results = runner.invoke(coverage.main, ["--help"])
    assert results.exit_code == 0

    results = runner.invoke(coverage.main, ["--version"])
    assert results.exit_code == 0


@pytest.mark.slow
def test_main(tmpdir):
    directory_run = tmpdir.mkdir("report")
    runner = CliRunner()

    results = runner.invoke(
        coverage.main,
        [
            "-i",
            bedfile,
            "-o",
            "--output-directory",
            str(directory_run),
            "--no-multiqc",
            "--window-median",
            "3001",
            "-r",
            fastafile,
        ],
    )
    assert os.path.exists(str(directory_run / "JB409847/sequana_summary_coverage.json"))
    assert results.exit_code == 0


def __test_main_downloads():
    runner = CliRunner()

    with runner.isolated_filesystem():

        results = runner.invoke(
            coverage.main,
            [
                "--download-reference",
                "JB409847",
            ],
        )
        assert results.exit_code == 0
        assert os.path.exists("JB409847.fa")

        results = runner.invoke(coverage.main, ["--download-genbank", "JB409847"])
        assert results.exit_code == 0
        assert os.path.exists("JB409847.gbk")


@pytest.mark.slow
def test_multiqc_report_and_annotation(tmpdir):
    directory_run = tmpdir.mkdir("report")
    runner = CliRunner()

    results = runner.invoke(
        coverage.main,
        [
            "-i",
            bedfile,
            "-o",
            "--output-directory",
            str(directory_run),
            "--window-median",
            "3001",
            "--annotation-file",
            genbank,
            "-r",
            fastafile,
            "-H",
            4,
            "-L",
            -4,
        ],
    )
    assert os.path.exists(str(directory_run / "JB409847/sequana_summary_coverage.json"))
    assert os.path.exists(str(directory_run / "multiqc_report.html"))
    assert results.exit_code == 0


def test_no_stray_directories_in_working_directory(tmp_path, monkeypatch):
    """Everything must be created in --output-directory, not in the current directory."""
    workdir = tmp_path / "work"
    workdir.mkdir()
    outdir = tmp_path / "elsewhere" / "report"
    monkeypatch.chdir(workdir)

    results = CliRunner().invoke(
        coverage.main,
        ["-i", bedfile, "--output-directory", str(outdir), "--no-multiqc", "--window-median", "3001"],
    )
    assert results.exit_code == 0, results.output
    assert (outdir / "JB409847" / "subplots").is_dir()
    assert list(workdir.iterdir()) == []


def test_reference_without_matching_names_does_not_crash(tmp_path, caplog):
    """A reference that is not a FASTA of the right chromosomes (e.g. a GFF file) only warns."""
    wrong = tmp_path / "other.fa"
    wrong.write_text(">another_name\nACGTACGTACGT\n")
    results = CliRunner().invoke(
        coverage.main,
        [
            "-i",
            bedfile,
            "--output-directory",
            str(tmp_path / "report"),
            "--no-multiqc",
            "--window-median",
            "3001",
            "--reference-file",
            str(wrong),
        ],
    )
    assert results.exit_code == 0, results.output
    assert (tmp_path / "report" / "JB409847" / "JB409847.cov.html").exists()


def test_viewer_option(tmp_path):
    """the viewer (default) adds an interactive page per chromosome, linked from the chromosome report."""
    import json
    import re

    outdir = tmp_path / "report"
    results = CliRunner().invoke(
        coverage.main,
        [
            "-i",
            bedfile,
            "--output-directory",
            str(outdir),
            "--no-multiqc",
            "--window-median",
            "3001",
            "--reference-file",
            fastafile,
            "--annotation-file",
            genbank,
        ],
    )
    assert results.exit_code == 0, results.output
    page = outdir / "JB409847" / "JB409847.viewer.html"
    assert page.exists()
    text = page.read_text()
    meta = json.loads(
        re.search(r'<script id="meta" type="application/json">(.*?)</script>', text, re.S)
        .group(1)
        .replace("<\\/", "</")
    )
    assert meta["name"] == "JB409847" and meta["has_gc"] and meta["thresholds"]["high"] == 4
    assert len(meta["ann"]["features"]) > 0 and len(meta["rois"]) > 0  # genes from the GenBank file
    assert "JB409847.viewer.html" in (outdir / "JB409847" / "JB409847.cov.html").read_text()

    # explicitly disabled: no viewer
    plain = tmp_path / "plain"
    results = CliRunner().invoke(
        coverage.main,
        ["-i", bedfile, "--output-directory", str(plain), "--no-multiqc", "--window-median", "3001", "--no-viewer"],
    )
    assert results.exit_code == 0, results.output
    assert not (plain / "JB409847" / "JB409847.viewer.html").exists()
    assert "viewer.html" not in (plain / "JB409847" / "JB409847.cov.html").read_text()


def test_engine_rust(tmp_path):
    """--engine rust finds the same ROIs as the python engine (skipped without the extension)."""
    pytest.importorskip("sequana_cov_rs")
    rois = {}
    for engine in ("python", "rust"):
        outdir = tmp_path / engine
        results = CliRunner().invoke(
            coverage.main,
            ["-i", bedfile, "--output-directory", str(outdir), "--no-multiqc", "--window-median", "3001"]
            + ["--engine", engine],
        )
        assert results.exit_code == 0, results.output
        rois[engine] = pd.read_csv(outdir / "JB409847" / "rois.csv", index_col=0)
    # same regions; the fitted model differs in the last digits only
    pd.testing.assert_frame_equal(rois["python"], rois["rust"], check_dtype=False, rtol=1e-6)


def test_engine_rust_missing_extension(tmp_path, monkeypatch):
    import sequana.bedtools_rs as rs

    monkeypatch.setattr(rs, "rust_available", lambda: False)
    results = CliRunner().invoke(
        coverage.main, ["-i", bedfile, "--output-directory", str(tmp_path / "x"), "--engine", "rust"]
    )
    assert results.exit_code != 0
    assert "sequana_cov_rs" in results.output
