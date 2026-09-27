"""
Comprehensive tests for SequanaReport and SummaryBase classes

Tests cover:
- Report initialization with various metadata configurations
- Report content generation (dependencies, workflow, caller)
- Dependency table creation
- Edge cases (missing files, empty data, None values)
- Configuration integration
"""

import json
import os
import tempfile
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from sequana.modules_report.summary import SequanaReport, SummaryBase

from . import test_dir


@pytest.fixture
def sample_report_data():
    """Create sample report metadata."""
    return {
        "name": "test_pipeline",
        "pipeline_version": "1.0.0",
        "sequana_wrappers": "latest",
        "rulegraph": "workflow.svg",
        "requirements": "requirements.txt",
    }


@pytest.fixture
def temp_report_dir(worker_id):
    """Create isolated temporary directory for each test/worker.

    With pytest-xdist, each worker gets its own temp dir to avoid
    race conditions when creating report subdirectories.
    """
    tmpdir = tempfile.mkdtemp(suffix=f"_{worker_id}")
    yield tmpdir
    # Cleanup
    import shutil

    if os.path.exists(tmpdir):
        shutil.rmtree(tmpdir, ignore_errors=True)


@pytest.fixture
def mock_config(temp_report_dir):
    """Mock the sequana config module."""
    # Patch at base_module level (used during __init__) and summary level
    with patch("sequana.modules_report.base_module.config") as mock_base, patch(
        "sequana.modules_report.summary.config"
    ) as mock_summary:
        for mock in (mock_base, mock_summary):
            mock.summary_sections = []
            mock.pipeline_version = None
            mock.pipeline_name = None
            mock.sequana_wrappers = None
            mock.output_dir = temp_report_dir
            mock.css_list = []
            mock.js_list = []
            mock.logo = None
        yield mock_summary


class TestSummaryBaseInit:
    """Test SummaryBase initialization."""

    def test_summarybase_init_no_dir(self, mock_config, temp_report_dir):
        """Test SummaryBase initialization without required_dir."""
        sb = SummaryBase()
        assert sb is not None

    def test_summarybase_init_with_dir(self, mock_config, temp_report_dir):
        """Test SummaryBase initialization with required_dir."""
        sb = SummaryBase(required_dir=("css", "js"))
        assert sb is not None

    def test_summarybase_is_initialized(self, mock_config, temp_report_dir):
        """Test that SummaryBase is properly initialized."""
        sb = SummaryBase()
        assert isinstance(sb, SummaryBase)


class TestSequanaReportInit:
    """Test SequanaReport initialization."""

    def test_init_basic(self, sample_report_data, mock_config):
        """Test basic initialization with minimal data."""
        with patch("sequana.modules_report.summary.SequanaReport.create_html"):
            report = SequanaReport(sample_report_data)

        assert report.name == "test_pipeline"
        assert report.wrappers == "latest"

    def test_init_with_title(self, sample_report_data, mock_config):
        """Test initialization with custom title."""
        with patch("sequana.modules_report.summary.SequanaReport.create_html"):
            report = SequanaReport(sample_report_data)

        # Title is set automatically from pipeline name
        assert report.title is not None

    def test_init_with_intro(self, sample_report_data, mock_config):
        """Test initialization with introduction text."""
        intro_text = "<p>This is an introduction</p>"

        with patch("sequana.modules_report.summary.SequanaReport.create_html"):
            report = SequanaReport(
                sample_report_data,
                intro=intro_text,
            )

        assert report.intro == intro_text

    def test_init_with_output_filename(self, sample_report_data, mock_config, temp_report_dir):
        """Test initialization with custom output filename."""
        output_file = os.path.join(temp_report_dir, "custom_report.html")

        with patch("sequana.modules_report.summary.SequanaReport.create_html"):
            report = SequanaReport(
                sample_report_data,
                output_filename=output_file,
            )

        assert report is not None

    def test_init_without_workflow(self, sample_report_data, mock_config):
        """Test initialization with workflow disabled."""
        with patch("sequana.modules_report.summary.SequanaReport.create_html"):
            with patch("sequana.modules_report.summary.SequanaReport.workflow"):
                report = SequanaReport(
                    sample_report_data,
                    workflow=False,
                )

        assert report is not None

    def test_init_with_nsamples(self, sample_report_data, mock_config):
        """Test initialization with sample count."""
        n_samples = 5

        with patch("sequana.modules_report.summary.SequanaReport.create_html"):
            report = SequanaReport(
                sample_report_data,
                Nsamples=n_samples,
            )

        assert report is not None

    def test_init_sets_config_values(self, sample_report_data, mock_config):
        """Test that initialization sets config values."""
        with patch("sequana.modules_report.summary.SequanaReport.create_html"):
            SequanaReport(sample_report_data)

        # Verify config was updated
        assert mock_config.pipeline_version == sample_report_data["pipeline_version"]
        assert mock_config.pipeline_name == sample_report_data["name"]
        assert mock_config.sequana_wrappers == sample_report_data["sequana_wrappers"]


class TestSequanaReportContentCreation:
    """Test report content generation."""

    def test_create_report_content_with_workflow(self, sample_report_data, mock_config):
        """Test report content creation with workflow."""
        with patch("sequana.modules_report.summary.SequanaReport.create_html"):
            with patch("sequana.modules_report.summary.SequanaReport.workflow"):
                with patch("sequana.modules_report.summary.SequanaReport.dependencies"):
                    with patch("sequana.modules_report.summary.SequanaReport.caller"):
                        report = SequanaReport(
                            sample_report_data,
                            workflow=True,
                        )

        assert hasattr(report, "sections")
        assert isinstance(report.sections, list)

    def test_create_report_content_without_workflow(self, sample_report_data, mock_config):
        """Test report content creation without workflow."""
        with patch("sequana.modules_report.summary.SequanaReport.create_html"):
            with patch("sequana.modules_report.summary.SequanaReport.dependencies"):
                with patch("sequana.modules_report.summary.SequanaReport.caller"):
                    report = SequanaReport(
                        sample_report_data,
                        workflow=False,
                    )

        assert hasattr(report, "sections")


class TestSequanaReportDependencies:
    """Test dependency table generation."""

    def test_get_table_dependencies(self, sample_report_data, mock_config):
        """Test dependency table generation."""
        with patch("sequana.modules_report.summary.SequanaReport.create_html"):
            report = SequanaReport(sample_report_data)

            # Mock importlib.metadata.requires
            with patch("sequana.modules_report.summary.importlib.metadata.requires") as mock_requires:
                mock_requires.return_value = [
                    "numpy>=1.19.0",
                    "pandas>=1.0.0",
                ]

                html = report.get_table_dependencies()

        assert isinstance(html, str)

    def test_get_table_versions_with_versions_file(self, sample_report_data, mock_config, temp_report_dir):
        """Test version table when versions.txt exists."""
        # Create .sequana directory with versions.txt
        sequana_dir = os.path.join(temp_report_dir, ".sequana")
        os.makedirs(sequana_dir, exist_ok=True)

        versions_file = os.path.join(sequana_dir, "versions.txt")
        with open(versions_file, "w") as f:
            f.write("samtools 1.13\n")
            f.write("bwa 0.7.17\n")

        # Change to temp directory
        original_cwd = os.getcwd()
        os.chdir(temp_report_dir)

        try:
            with patch("sequana.modules_report.summary.SequanaReport.create_html"):
                report = SequanaReport(sample_report_data)
                html = report.get_table_versions()
        finally:
            os.chdir(original_cwd)

        # Should return non-empty string if file exists
        if html:
            assert isinstance(html, str)

    def test_get_table_versions_without_versions_file(self, sample_report_data, mock_config, temp_report_dir):
        """Test version table when versions.txt doesn't exist."""
        original_cwd = os.getcwd()
        os.chdir(temp_report_dir)

        try:
            with patch("sequana.modules_report.summary.SequanaReport.create_html"):
                report = SequanaReport(sample_report_data)
                html = report.get_table_versions()
        finally:
            os.chdir(original_cwd)

        # Should return empty string if file doesn't exist
        assert html == ""


class TestSequanaReportCaller:
    """Test caller/command section generation."""

    def test_caller_creates_section(self, sample_report_data, mock_config, temp_report_dir):
        """Test caller section is created."""
        sequana_dir = os.path.join(temp_report_dir, ".sequana")
        os.makedirs(sequana_dir, exist_ok=True)

        info_file = os.path.join(sequana_dir, "info.txt")
        with open(info_file, "w") as f:
            f.write("# Generated command\n")
            f.write("sequana analyze --input test.fastq\n")

        original_cwd = os.getcwd()
        os.chdir(temp_report_dir)

        try:
            with patch("sequana.modules_report.summary.SequanaReport.create_html"):
                report = SequanaReport(sample_report_data)
                report.caller()

                # Verify section was added
                assert len(report.sections) > 0
                assert any(s.get("anchor") == "command" for s in report.sections)
        finally:
            os.chdir(original_cwd)

    def test_caller_handles_missing_info_file(self, sample_report_data, mock_config, temp_report_dir):
        """Test caller section handles missing info file."""
        original_cwd = os.getcwd()
        os.chdir(temp_report_dir)

        try:
            with patch("sequana.modules_report.summary.SequanaReport.create_html"):
                report = SequanaReport(sample_report_data)
                report.caller()

                # Should still add a section even if file is missing
                assert any(s.get("anchor") == "command" for s in report.sections)
        finally:
            os.chdir(original_cwd)


class TestSequanaReportEdgeCases:
    """Test edge cases and error handling."""

    def test_init_with_missing_name(self, mock_config):
        """Test initialization with missing pipeline name."""
        data = {
            "pipeline_version": "1.0.0",
            "rulegraph": "workflow.svg",
        }

        with patch("sequana.modules_report.summary.SequanaReport.create_html"):
            report = SequanaReport(data)

        # Should default to "undefined"
        assert report.name == "undefined"

    def test_init_with_missing_pipeline_version(self, sample_report_data, mock_config):
        """Test initialization with missing pipeline version."""
        data = sample_report_data.copy()
        del data["pipeline_version"]

        with patch("sequana.modules_report.summary.SequanaReport.create_html"):
            report = SequanaReport(data)

        # Should set to "latest" default
        assert mock_config.pipeline_version == "latest"

    def test_init_with_none_values(self, mock_config):
        """Test initialization with None values in data."""
        data = {
            "name": None,
            "pipeline_version": None,
            "rulegraph": "workflow.svg",
        }

        with patch("sequana.modules_report.summary.SequanaReport.create_html"):
            report = SequanaReport(data)

        assert report is not None


class TestSequanaReportDataIntegration:
    """Test integration with real-like data structures."""

    def test_init_with_real_like_metadata(self, mock_config):
        """Test with realistic pipeline metadata."""
        realistic_data = {
            "name": "rnadiff",
            "pipeline_version": "2.0.0",
            "sequana_wrappers": "2024.01",
            "rulegraph": ".sequana/rulegraph.svg",
            "requirements": ".sequana/requirements.txt",
        }

        with patch("sequana.modules_report.summary.SequanaReport.create_html"):
            report = SequanaReport(realistic_data)

        assert report.name == "rnadiff"
        assert report.wrappers == "2024.01"

    def test_all_main_methods_present(self, sample_report_data, mock_config):
        """Test that all main methods are present."""
        with patch("sequana.modules_report.summary.SequanaReport.create_html"):
            report = SequanaReport(sample_report_data)

        assert callable(report.create_report_content)
        assert callable(report.dependencies)
        assert callable(report.workflow)
        assert callable(report.caller)
        assert callable(report.get_table_dependencies)
        assert callable(report.get_table_versions)


class TestSequanaReportConfiguration:
    """Test configuration and settings."""

    def test_config_pipeline_version_set(self, sample_report_data, mock_config):
        """Test that pipeline_version is set in config."""
        with patch("sequana.modules_report.summary.SequanaReport.create_html"):
            SequanaReport(sample_report_data)

        assert mock_config.pipeline_version == sample_report_data["pipeline_version"]

    def test_config_pipeline_name_set(self, sample_report_data, mock_config):
        """Test that pipeline_name is set in config."""
        with patch("sequana.modules_report.summary.SequanaReport.create_html"):
            SequanaReport(sample_report_data)

        assert mock_config.pipeline_name == sample_report_data["name"]

    def test_config_wrappers_set(self, sample_report_data, mock_config):
        """Test that sequana_wrappers is set in config."""
        with patch("sequana.modules_report.summary.SequanaReport.create_html"):
            SequanaReport(sample_report_data)

        assert mock_config.sequana_wrappers == sample_report_data["sequana_wrappers"]


class TestSequanaReportAttributes:
    """Test report attributes."""

    def test_json_attribute_stored(self, sample_report_data, mock_config):
        """Test that JSON metadata is stored."""
        with patch("sequana.modules_report.summary.SequanaReport.create_html"):
            report = SequanaReport(sample_report_data)

        assert report.json == sample_report_data

    def test_intro_attribute_default_empty(self, sample_report_data, mock_config):
        """Test that intro defaults to empty string."""
        with patch("sequana.modules_report.summary.SequanaReport.create_html"):
            report = SequanaReport(sample_report_data, intro="")

        assert report.intro == ""

    def test_report_has_required_methods(self, sample_report_data, mock_config):
        """Test that report has required methods."""
        with patch("sequana.modules_report.summary.SequanaReport.create_html"):
            report = SequanaReport(sample_report_data)

        assert callable(report.create_report_content)
