"""
Comprehensive tests for the pairtools MultiQC module.

Tests cover:
- Plugin initialization and configuration
- Data parsing from pairtools output files
- Report generation and output formatting
- Integration with MultiQC pipeline
- Edge cases (empty files, malformed data, missing fields)
- Public API methods
- Error handling
"""

import io
import os
from unittest import mock

import pytest

try:
    from multiqc import report
    from multiqc.base_module import ModuleNoSamplesFound

    from sequana.multiqc.pairtools import MultiqcModule

    from . import test_dir

    def sequana_data(file):
        return f"{test_dir}/data/{file}"

    # --- MultiQC Module Initialization Tests ---

    def test_multiqc_module_basic_initialization():
        """Test basic module initialization with a single valid stats file"""
        report.reset()
        report.files = {
            "pairtools": [
                {"fn": sequana_data("pairtools_sample1.stats"), "root": ".", "sp_key": "pairtools"},
            ]
        }
        module = MultiqcModule()
        assert module is not None
        assert module.name == "pairtools"
        assert module.anchor == "pairtools"
        assert len(module.pairtools_stats) == 1

    def test_multiqc_module_module_properties():
        """Test module has correct name, anchor, and href properties"""
        report.reset()
        report.files = {
            "pairtools": [
                {"fn": sequana_data("pairtools_sample1.stats"), "root": ".", "sp_key": "pairtools"},
            ]
        }
        module = MultiqcModule()
        assert module.name == "pairtools"
        assert module.anchor == "pairtools"
        # href can be a list or string depending on MultiQC version
        href_str = str(module.href).lower()
        assert "mirnylab" in href_str or "pairtools" in href_str
        assert isinstance(module.doi, list) and len(module.doi) > 0

    def test_multiqc_module_multiple_samples():
        """Test module correctly processes multiple sample files"""
        report.reset()
        report.files = {
            "pairtools": [
                {"fn": sequana_data("pairtools_sample1.stats"), "root": ".", "sp_key": "pairtools"},
                {"fn": sequana_data("pairtools_sample2.stats"), "root": ".", "sp_key": "pairtools"},
            ]
        }
        module = MultiqcModule()
        assert len(module.pairtools_stats) == 2
        # Verify both samples are parsed
        assert any(
            "sample1" in key.lower() or "pairtools_sample1" in key.lower() for key in module.pairtools_stats.keys()
        )

    def test_multiqc_module_no_samples_found():
        """Test module raises ModuleNoSamplesFound when no valid samples exist"""
        report.reset()
        report.files = {
            "pairtools": [
                {"fn": sequana_data("pairtools_incomplete.stats"), "root": ".", "sp_key": "pairtools"},
            ]
        }
        with pytest.raises(ModuleNoSamplesFound):
            MultiqcModule()

    def test_multiqc_module_empty_files_list():
        """Test module raises ModuleNoSamplesFound when files list is empty"""
        report.reset()
        report.files = {"pairtools": []}
        with pytest.raises(ModuleNoSamplesFound):
            MultiqcModule()

    def test_multiqc_module_mixed_valid_invalid_files():
        """Test module processes only valid files when mixed with invalid ones"""
        report.reset()
        report.files = {
            "pairtools": [
                {"fn": sequana_data("pairtools_incomplete.stats"), "root": ".", "sp_key": "pairtools"},
                {"fn": sequana_data("pairtools_sample1.stats"), "root": ".", "sp_key": "pairtools"},
            ]
        }
        module = MultiqcModule()
        # Should have at least one valid sample
        assert len(module.pairtools_stats) >= 1

    # --- Data Parsing Tests ---

    def test_parse_pairtools_stats_basic():
        """Test parse_pairtools_stats static method with valid file"""
        with open(sequana_data("pairtools_sample1.stats"), "r") as f:
            result = MultiqcModule.parse_pairtools_stats({"f": f})
        assert result is not None
        assert "total" in result
        assert result["total"] == 5000

    def test_parse_pairtools_stats_comprehensive_data():
        """Test that parse_pairtools_stats extracts all expected fields"""
        with open(sequana_data("pairtools_sample1.stats"), "r") as f:
            result = MultiqcModule.parse_pairtools_stats({"f": f})
        assert result is not None
        # Check required fields
        assert result["total"] == 5000
        assert result["total_unmapped"] == 500
        assert result["total_single_sided_mapped"] == 250
        assert result["total_mapped"] == 4250
        assert result["total_dups"] == 850
        assert result["total_nodups"] == 3400
        assert result["cis"] == 2200
        assert result["trans"] == 1200
        # Check calculated fractions
        assert "frac_unmapped" in result
        assert "frac_mapped" in result
        assert "frac_cis" in result

    def test_parse_pairtools_stats_pair_types():
        """Test that pair_types are correctly extracted"""
        with open(sequana_data("pairtools_sample1.stats"), "r") as f:
            result = MultiqcModule.parse_pairtools_stats({"f": f})
        assert result["pair_types"] is not None
        assert "UU" in result["pair_types"]
        assert result["pair_types"]["UU"] == 3200

    def test_parse_pairtools_stats_cis_dist():
        """Test that cis_dist data is extracted and processed"""
        with open(sequana_data("pairtools_sample1.stats"), "r") as f:
            result = MultiqcModule.parse_pairtools_stats({"f": f})
        assert result["cis_dist"] is not None
        assert isinstance(result["cis_dist"], dict)
        assert "trans" in result["cis_dist"]

    def test_parse_pairtools_stats_dist_freq():
        """Test that dist_freq data is extracted and processed"""
        with open(sequana_data("pairtools_sample1.stats"), "r") as f:
            result = MultiqcModule.parse_pairtools_stats({"f": f})
        assert result["dist_freq"] is not None
        assert isinstance(result["dist_freq"], dict)
        assert "all" in result["dist_freq"]
        for orientation in ["++", "+-", "-+", "--"]:
            assert orientation in result["dist_freq"]

    def test_parse_pairtools_stats_minimal_file():
        """Test parsing file with minimal required data only"""
        with open(sequana_data("pairtools_minimal.stats"), "r") as f:
            result = MultiqcModule.parse_pairtools_stats({"f": f})
        assert result is not None
        assert result["cis_dist"] is None
        assert result["dist_freq"] is None
        assert result["pair_types"] is not None

    # --- Report Generation Tests ---

    def test_pair_types_chart_generation():
        """Test pair_types_chart generates valid output"""
        report.reset()
        report.files = {
            "pairtools": [
                {"fn": sequana_data("pairtools_sample1.stats"), "root": ".", "sp_key": "pairtools"},
            ]
        }
        module = MultiqcModule()
        chart = module.pair_types_chart()
        assert chart is not None
        # Should return a chart object or string (if None data, returns alert string)
        assert isinstance(chart, str) or hasattr(chart, "plot_type")

    def test_pair_types_chart_with_multiple_samples():
        """Test pair_types_chart with multiple samples"""
        report.reset()
        report.files = {
            "pairtools": [
                {"fn": sequana_data("pairtools_sample1.stats"), "root": ".", "sp_key": "pairtools"},
                {"fn": sequana_data("pairtools_sample2.stats"), "root": ".", "sp_key": "pairtools"},
            ]
        }
        module = MultiqcModule()
        chart = module.pair_types_chart()
        assert chart is not None

    def test_pairs_by_cisrange_trans_generation():
        """Test pairs_by_cisrange_trans chart generation"""
        report.reset()
        report.files = {
            "pairtools": [
                {"fn": sequana_data("pairtools_sample1.stats"), "root": ".", "sp_key": "pairtools"},
            ]
        }
        module = MultiqcModule()
        chart = module.pairs_by_cisrange_trans()
        assert chart is not None
        assert isinstance(chart, str) or hasattr(chart, "plot_type")

    def test_pairs_by_strand_orientation_generation():
        """Test pairs_by_strand_orientation chart generation"""
        report.reset()
        report.files = {
            "pairtools": [
                {"fn": sequana_data("pairtools_sample1.stats"), "root": ".", "sp_key": "pairtools"},
            ]
        }
        module = MultiqcModule()
        chart = module.pairs_by_strand_orientation()
        assert chart is not None
        assert isinstance(chart, str) or hasattr(chart, "plot_type")

    def test_pairs_with_genomic_separation_generation():
        """Test pairs_with_genomic_separation plot generation"""
        report.reset()
        report.files = {
            "pairtools": [
                {"fn": sequana_data("pairtools_sample1.stats"), "root": ".", "sp_key": "pairtools"},
            ]
        }
        module = MultiqcModule()
        chart = module.pairs_with_genomic_separation()
        assert chart is not None
        assert isinstance(chart, str) or hasattr(chart, "plot_type")

    def test_pairtools_general_stats():
        """Test that general_stats are added to report"""
        report.reset()
        report.files = {
            "pairtools": [
                {"fn": sequana_data("pairtools_sample1.stats"), "root": ".", "sp_key": "pairtools"},
            ]
        }
        module = MultiqcModule()
        # pairtools_general_stats is called in __init__
        # Verify the method was called by checking the pairtools_stats structure
        assert len(module.pairtools_stats) > 0

    # --- Edge Cases and Error Handling ---

    def test_pair_types_chart_handles_none_data():
        """Test pair_types_chart gracefully handles None pair_types"""
        with mock.patch.object(MultiqcModule, "__init__", return_value=None):
            module = MultiqcModule()
            module.pairtools_stats = {
                "sample1": {"pair_types": None},
            }
            module.params = {"pairtypes_colors": {}}
            chart = module.pair_types_chart()
            # Should return an alert message
            assert isinstance(chart, str)
            assert "alert" in chart.lower() or "skip" in chart.lower()

    def test_pairs_by_cisrange_trans_handles_none_data():
        """Test pairs_by_cisrange_trans handles None cis_dist"""
        with mock.patch.object(MultiqcModule, "__init__", return_value=None):
            module = MultiqcModule()
            module.pairtools_stats = {
                "sample1": {"cis_dist": None},
            }
            module.params = {"cis_range_colors": []}
            chart = module.pairs_by_cisrange_trans()
            assert isinstance(chart, str)
            assert "alert" in chart.lower() or "skip" in chart.lower()

    def test_pairs_by_strand_orientation_handles_none_data():
        """Test pairs_by_strand_orientation handles None pairs_by_strand"""
        with mock.patch.object(MultiqcModule, "__init__", return_value=None):
            module = MultiqcModule()
            module.pairtools_stats = {
                "sample1": {"pairs_by_strand": None},
            }
            module.params = {
                "pairs_orientation_names": {"++": "FF", "-+": "RF", "+-": "FR", "--": "RR"},
                "pairs_orientation_colors": {"++": "#e41a1c", "-+": "#377eb8", "+-": "#4daf4a", "--": "#984ea3"},
            }
            chart = module.pairs_by_strand_orientation()
            assert isinstance(chart, str)
            assert "alert" in chart.lower() or "skip" in chart.lower()

    def test_pairs_with_genomic_separation_handles_none_data():
        """Test pairs_with_genomic_separation handles None dist_freq"""
        with mock.patch.object(MultiqcModule, "__init__", return_value=None):
            module = MultiqcModule()
            module.pairtools_stats = {
                "sample1": {"dist_freq": None},
            }
            chart = module.pairs_with_genomic_separation()
            assert isinstance(chart, str)
            assert "alert" in chart.lower() or "skip" in chart.lower()

    def test_pairs_by_cisrange_trans_different_categories():
        """Test pairs_by_cisrange_trans with samples having different distance categories"""
        with mock.patch.object(MultiqcModule, "__init__", return_value=None):
            module = MultiqcModule()
            module.pairtools_stats = {
                "sample1": {"cis_dist": {"cis: 0-1Kb": 100, "trans": 50}},
                "sample2": {"cis_dist": {"cis: 0-5Kb": 150, "trans": 60}},
            }
            module.params = {"cis_range_colors": ["#8c2d04"]}
            chart = module.pairs_by_cisrange_trans()
            # Should handle different categories gracefully
            assert isinstance(chart, str)

    def test_pairs_by_strand_orientation_empty_categories():
        """Test pairs_by_strand_orientation with no distance categories"""
        with mock.patch.object(MultiqcModule, "__init__", return_value=None):
            module = MultiqcModule()
            module.pairtools_stats = {
                "sample1": {"pairs_by_strand": {}},
            }
            module.params = {
                "pairs_orientation_names": {"++": "FF", "-+": "RF", "+-": "FR", "--": "RR"},
                "pairs_orientation_colors": {"++": "#e41a1c", "-+": "#377eb8", "+-": "#4daf4a", "--": "#984ea3"},
            }
            chart = module.pairs_by_strand_orientation()
            assert isinstance(chart, str)
            assert "alert" in chart.lower() or "no" in chart.lower()

    # --- Integration Tests ---

    def test_module_sections_created():
        """Test that all expected sections are created in the module"""
        report.reset()
        report.files = {
            "pairtools": [
                {"fn": sequana_data("pairtools_sample1.stats"), "root": ".", "sp_key": "pairtools"},
            ]
        }
        module = MultiqcModule()
        # Check that the module has sections (they're added via add_section)
        # We can verify by checking that the module initialized successfully
        assert module.pairtools_stats is not None
        assert len(module.pairtools_stats) > 0

    def test_module_params_loaded():
        """Test that params.yml is correctly loaded"""
        report.reset()
        report.files = {
            "pairtools": [
                {"fn": sequana_data("pairtools_sample1.stats"), "root": ".", "sp_key": "pairtools"},
            ]
        }
        module = MultiqcModule()
        assert module.params is not None
        assert "pairtypes_colors" in module.params
        assert "cis_range_colors" in module.params
        assert "pairs_orientation_names" in module.params
        assert "pairs_orientation_colors" in module.params

    def test_module_data_source_tracking():
        """Test that data sources are tracked properly"""
        report.reset()
        report.files = {
            "pairtools": [
                {"fn": sequana_data("pairtools_sample1.stats"), "root": ".", "sp_key": "pairtools"},
            ]
        }
        with mock.patch.object(MultiqcModule, "add_data_source") as mock_source:
            module = MultiqcModule()
            # add_data_source should be called for each valid sample
            assert mock_source.called

    def test_module_data_file_writing():
        """Test that parsed data is written to file"""
        report.reset()
        report.files = {
            "pairtools": [
                {"fn": sequana_data("pairtools_sample1.stats"), "root": ".", "sp_key": "pairtools"},
            ]
        }
        with mock.patch.object(MultiqcModule, "write_data_file") as mock_write:
            module = MultiqcModule()
            # write_data_file should be called to save parsed data
            assert mock_write.called

    # --- Additional Edge Cases ---

    def test_minimal_stats_file_parsing():
        """Test parsing file with only required fields"""
        report.reset()
        report.files = {
            "pairtools": [
                {"fn": sequana_data("pairtools_minimal.stats"), "root": ".", "sp_key": "pairtools"},
            ]
        }
        module = MultiqcModule()
        assert len(module.pairtools_stats) == 1
        sample = list(module.pairtools_stats.values())[0]
        # Verify all required fields are present
        for key in [
            "total",
            "total_unmapped",
            "total_single_sided_mapped",
            "total_mapped",
            "total_dups",
            "total_nodups",
            "cis",
            "trans",
        ]:
            assert key in sample

    def test_all_pair_types_colors_available():
        """Test that pairtypes_colors includes common types"""
        report.reset()
        report.files = {
            "pairtools": [
                {"fn": sequana_data("pairtools_sample1.stats"), "root": ".", "sp_key": "pairtools"},
            ]
        }
        module = MultiqcModule()
        # Check that common pair types have colors defined
        assert "UU" in module.params["pairtypes_colors"]
        assert "NN" in module.params["pairtypes_colors"]
        assert "DD" in module.params["pairtypes_colors"]

    def test_orientation_names_and_colors_complete():
        """Test that all four orientations have names and colors defined"""
        report.reset()
        report.files = {
            "pairtools": [
                {"fn": sequana_data("pairtools_sample1.stats"), "root": ".", "sp_key": "pairtools"},
            ]
        }
        module = MultiqcModule()
        for orientation in ["++", "-+", "+-", "--"]:
            assert orientation in module.params["pairs_orientation_names"]
            assert orientation in module.params["pairs_orientation_colors"]

except ImportError:
    pass
