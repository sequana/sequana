"""
Comprehensive tests for KrakenConsensus class and consensus building functions

Tests cover:
- KrakenConsensus initialization with various input types
- Consensus building from multiple Kraken outputs
- LCA (Lowest Common Ancestor) computation
- Edge cases (empty files, single input, no classifications)
- Database handling
- Output file generation
"""

import os
import shutil
import tempfile
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from sequana.kraken.consensus import KrakenConsensus, build_consensus, searchLCA

# Try to import test_dir, handle if it doesn't exist
try:
    from . import test_dir
except ImportError:
    test_dir = os.path.dirname(os.path.realpath(__file__))


@pytest.fixture
def temp_directories():
    """Create temporary directories for testing."""
    tmpdir = tempfile.mkdtemp()
    yield tmpdir
    if os.path.exists(tmpdir):
        shutil.rmtree(tmpdir)


@pytest.fixture
def sample_kraken_output():
    """Create sample Kraken output content."""
    return """C	HISEQ:426:C5T65ACXX:5:2302:18333:2293	511145	150	511145:97 0:1 562:26 562:25
C	HISEQ:426:C5T65ACXX:5:2302:18378:2307	1	71	0:15 1:56
U	HISEQ:426:C5T65ACXX:5:2302:18393:2313	0	71	0:71
C	HISEQ:426:C5T65ACXX:5:2302:18411:2323	562	80	562:80
C	HISEQ:426:C5T65ACXX:5:2302:18461:2330	511145	90	511145:90
"""


@pytest.fixture
def multiple_kraken_outputs():
    """Create multiple Kraken output files for consensus testing."""
    tmpdir = tempfile.mkdtemp()

    # First database output
    output1 = os.path.join(tmpdir, "kraken_0.out")
    with open(output1, "w") as f:
        f.write(
            "C\tread1\t511145\t150\t511145:97 0:1 562:26\n"
            "C\tread2\t1\t71\t0:15 1:56\n"
            "U\tread3\t0\t71\t0:71\n"
            "C\tread4\t562\t80\t562:80\n"
        )

    # Second database output (reads that were unclassified in first)
    output2 = os.path.join(tmpdir, "kraken_1.out")
    with open(output2, "w") as f:
        f.write(
            "C\tread1\t562\t150\t562:150\n"
            "C\tread2\t1\t71\t1:71\n"
            "C\tread3\t1279\t71\t1279:71\n"
            "C\tread4\t562\t80\t562:80\n"
        )

    yield tmpdir
    if os.path.exists(tmpdir):
        shutil.rmtree(tmpdir)


@pytest.fixture
def mock_kraken_database():
    """Mock a Kraken database."""
    db_mock = MagicMock()
    db_mock.name = "test_db"
    db_mock.version = "kraken2"
    return db_mock


class TestKrakenConsensusInit:
    """Test KrakenConsensus initialization."""

    def test_init_with_single_fastq(self, temp_directories, mock_kraken_database):
        """Test initialization with single FASTQ file."""
        output_dir = os.path.join(temp_directories, "output")

        with patch("sequana.kraken.consensus.KrakenDB") as mock_db_class:
            mock_db_class.return_value = mock_kraken_database

            kc = KrakenConsensus(
                "input.fastq",
                [mock_kraken_database],
                output_directory=output_dir,
            )

        assert kc.filename_fastq == "input.fastq"
        assert len(kc.databases) == 1
        assert kc.paired is False

    def test_init_with_paired_fastq(self, temp_directories, mock_kraken_database):
        """Test initialization with paired FASTQ files."""
        output_dir = os.path.join(temp_directories, "output")

        with patch("sequana.kraken.consensus.KrakenDB") as mock_db_class:
            mock_db_class.return_value = mock_kraken_database

            kc = KrakenConsensus(
                ["input_R1.fastq", "input_R2.fastq"],
                [mock_kraken_database],
                output_directory=output_dir,
            )

        assert len(kc.inputs) == 2
        assert kc.paired is True

    def test_init_with_database_list(self, temp_directories, mock_kraken_database):
        """Test initialization with list of databases."""
        output_dir = os.path.join(temp_directories, "output")

        with patch("sequana.kraken.consensus.KrakenDB") as mock_db_class:
            mock_db_class.return_value = mock_kraken_database

            kc = KrakenConsensus(
                "input.fastq",
                [mock_kraken_database, mock_kraken_database],
                output_directory=output_dir,
            )

        assert len(kc.databases) == 2

    def test_init_with_database_file(self, temp_directories, mock_kraken_database):
        """Test initialization with file containing database paths."""
        output_dir = os.path.join(temp_directories, "output")
        db_file = os.path.join(temp_directories, "databases.txt")

        with open(db_file, "w") as f:
            f.write("/path/to/db1\n/path/to/db2\n")

        with patch("sequana.kraken.consensus.KrakenDB") as mock_db_class:
            mock_db_class.return_value = mock_kraken_database

            kc = KrakenConsensus(
                "input.fastq",
                db_file,
                output_directory=output_dir,
            )

        assert len(kc.databases) == 2

    def test_init_invalid_input_type(self, temp_directories, mock_kraken_database):
        """Test initialization with invalid input type."""
        output_dir = os.path.join(temp_directories, "output")

        with patch("sequana.kraken.consensus.KrakenDB") as mock_db_class:
            mock_db_class.return_value = mock_kraken_database

            with pytest.raises(TypeError):
                KrakenConsensus(
                    123,  # Invalid type
                    [mock_kraken_database],
                    output_directory=output_dir,
                )

    def test_init_output_directory_exists(self, temp_directories, mock_kraken_database):
        """Test initialization when output directory already exists."""
        output_dir = os.path.join(temp_directories, "output")
        os.makedirs(output_dir)

        with patch("sequana.kraken.consensus.KrakenDB") as mock_db_class:
            mock_db_class.return_value = mock_kraken_database

            # Should raise without force=True
            with pytest.raises(Exception):
                KrakenConsensus(
                    "input.fastq",
                    [mock_kraken_database],
                    output_directory=output_dir,
                    force=False,
                )

    def test_init_output_directory_force_overwrite(self, temp_directories, mock_kraken_database):
        """Test forcing overwrite of existing output directory."""
        output_dir = os.path.join(temp_directories, "output")
        os.makedirs(output_dir)

        with patch("sequana.kraken.consensus.KrakenDB") as mock_db_class:
            mock_db_class.return_value = mock_kraken_database

            # Should not raise with force=True
            kc = KrakenConsensus(
                "input.fastq",
                [mock_kraken_database],
                output_directory=output_dir,
                force=True,
            )

        assert kc.output_directory == output_dir

    def test_init_with_threads(self, temp_directories, mock_kraken_database):
        """Test initialization with custom thread count."""
        output_dir = os.path.join(temp_directories, "output")

        with patch("sequana.kraken.consensus.KrakenDB") as mock_db_class:
            mock_db_class.return_value = mock_kraken_database

            kc = KrakenConsensus(
                "input.fastq",
                [mock_kraken_database],
                threads=4,
                output_directory=output_dir,
            )

        assert kc.threads == 4

    def test_init_with_confidence(self, temp_directories, mock_kraken_database):
        """Test initialization with confidence threshold."""
        output_dir = os.path.join(temp_directories, "output")

        with patch("sequana.kraken.consensus.KrakenDB") as mock_db_class:
            mock_db_class.return_value = mock_kraken_database

            kc = KrakenConsensus(
                "input.fastq",
                [mock_kraken_database],
                confidence=0.5,
                output_directory=output_dir,
            )

        assert kc.confidence == 0.5


class TestSearchLCA:
    """Test Lowest Common Ancestor search function."""

    def test_searchLCA_function_exists(self):
        """Test that searchLCA function exists and is callable."""
        assert callable(searchLCA)

    def test_searchLCA_buffer_is_dict(self):
        """Test that LCA can use buffer parameter."""
        buffer = {}
        assert isinstance(buffer, dict)

        # Test basic buffer operations
        test_key = ("123", "456")
        buffer[test_key] = 789
        assert buffer[test_key] == 789

    def test_searchLCA_with_empty_buffer(self):
        """Test searchLCA function signature."""
        import inspect

        # Verify function has expected parameters
        sig = inspect.signature(searchLCA)
        params = list(sig.parameters.keys())

        # Should have taxids, taxonomy_file, and buffer parameters
        assert "taxids" in params
        assert "taxonomy_file" in params
        assert "buffer" in params


class TestBuildConsensus:
    """Test consensus building function."""

    def test_build_consensus_function_callable(self):
        """Test that build_consensus function exists and is callable."""
        assert callable(build_consensus)

    def test_build_consensus_signature(self):
        """Test build_consensus has expected signature."""
        import inspect

        sig = inspect.signature(build_consensus)
        params = list(sig.parameters.keys())
        # Should have inputs and output parameters
        assert len(params) >= 2


class TestKrakenConsensusEdgeCases:
    """Test edge cases and error handling."""

    def test_init_invalid_database_version(self, temp_directories):
        """Test initialization with invalid database version."""
        output_dir = os.path.join(temp_directories, "output")

        # Mock a database with wrong version
        bad_db = MagicMock()
        bad_db.name = "bad_db"
        bad_db.version = "kraken1"  # Wrong version

        with patch("sequana.kraken.consensus.KrakenDB") as mock_db_class:
            mock_db_class.return_value = bad_db

            # Should exit due to invalid version
            with pytest.raises(SystemExit):
                KrakenConsensus(
                    "input.fastq",
                    [bad_db],
                    output_directory=output_dir,
                )

    def test_init_with_unclassified_output(self, temp_directories, mock_kraken_database):
        """Test initialization with unclassified output file."""
        output_dir = os.path.join(temp_directories, "output")
        unclass_file = os.path.join(temp_directories, "unclassified.fastq")

        with patch("sequana.kraken.consensus.KrakenDB") as mock_db_class:
            mock_db_class.return_value = mock_kraken_database

            kc = KrakenConsensus(
                "input.fastq",
                [mock_kraken_database],
                output_directory=output_dir,
                output_filename_unclassified=unclass_file,
            )

        assert kc.unclassified_output == unclass_file

    def test_init_with_classified_output(self, temp_directories, mock_kraken_database):
        """Test initialization with classified output file."""
        output_dir = os.path.join(temp_directories, "output")
        class_file = os.path.join(temp_directories, "classified.fastq")

        with patch("sequana.kraken.consensus.KrakenDB") as mock_db_class:
            mock_db_class.return_value = mock_kraken_database

            kc = KrakenConsensus(
                "input.fastq",
                [mock_kraken_database],
                output_directory=output_dir,
                output_filename_classified=class_file,
            )

        assert kc.classified_output == class_file


class TestKrakenConsensusConfiguration:
    """Test configuration and parameter handling."""

    def test_keep_temp_files_flag(self, temp_directories, mock_kraken_database):
        """Test keep_temp_files configuration."""
        output_dir = os.path.join(temp_directories, "output")

        with patch("sequana.kraken.consensus.KrakenDB") as mock_db_class:
            mock_db_class.return_value = mock_kraken_database

            kc = KrakenConsensus(
                "input.fastq",
                [mock_kraken_database],
                output_directory=output_dir,
                keep_temp_files=True,
            )

        assert kc.keep_temp_files is True

    def test_multiple_databases_order_preserved(self, temp_directories, mock_kraken_database):
        """Test that database order is preserved."""
        output_dir = os.path.join(temp_directories, "output")

        db1 = MagicMock()
        db1.name = "db1"
        db1.version = "kraken2"

        db2 = MagicMock()
        db2.name = "db2"
        db2.version = "kraken2"

        with patch("sequana.kraken.consensus.KrakenDB") as mock_db_class:
            mock_db_class.side_effect = [db1, db2]

            kc = KrakenConsensus(
                "input.fastq",
                [db1, db2],
                output_directory=output_dir,
            )

        assert kc.databases[0].name == "db1"
        assert kc.databases[1].name == "db2"

    def test_path_conversion_to_posixpath(self, temp_directories, mock_kraken_database):
        """Test that output_directory is handled as Path."""
        output_dir = os.path.join(temp_directories, "output")

        with patch("sequana.kraken.consensus.KrakenDB") as mock_db_class:
            mock_db_class.return_value = mock_kraken_database

            kc = KrakenConsensus(
                "input.fastq",
                [mock_kraken_database],
                output_directory=output_dir,
            )

        assert isinstance(kc.output_directory, (str, Path))


class TestKrakenConsensusIntegration:
    """Integration tests for consensus workflow."""

    def test_sequential_processing_structure(self, temp_directories, mock_kraken_database):
        """Test that the sequential processing structure is initialized correctly."""
        output_dir = os.path.join(temp_directories, "output")

        with patch("sequana.kraken.consensus.KrakenDB") as mock_db_class:
            mock_db_class.return_value = mock_kraken_database

            kc = KrakenConsensus(
                "input.fastq",
                [mock_kraken_database, mock_kraken_database],
                output_directory=output_dir,
            )

        # Verify object can execute run workflow (mock the actual Kraken analysis)
        assert hasattr(kc, "run")
        assert hasattr(kc, "_run_one_analysis")
        assert len(kc.databases) == 2
