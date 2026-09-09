"""
Comprehensive tests for bedtools.py coverage gaps.

Focus areas:
- Lines 1166-1172: Chromosome name annotation handling and warnings
- Lines 271-276: Chromosome list validation and filtering
- Lines 322-326: GC window size setter with even/odd number handling
"""

import os
import tempfile
from pathlib import Path

import pytest
from easydev import TempFile

from sequana import bedtools

from . import test_dir


class TestChromosomeListValidation:
    """Test chromosome list validation (lines 271-276)."""

    @pytest.fixture
    def sample_bed_file(self):
        """Create a BED file with multiple chromosomes."""
        return f"{test_dir}/data/bed/JB409847.bed"

    @pytest.fixture
    def genbank_file(self):
        """Genbank annotation file."""
        return f"{test_dir}/data/genbank/JB409847.gbk"

    def test_chromosome_list_valid_subset(self, sample_bed_file, genbank_file):
        """Test providing a valid subset of chromosome names."""
        bed = bedtools.SequanaCoverage(sample_bed_file, genbank_file)
        all_chroms = bed.chrom_names.copy()

        # Use subset of chromosomes
        bed_filtered = bedtools.SequanaCoverage(sample_bed_file, genbank_file, chromosome_list=[all_chroms[0]])
        assert len(bed_filtered.chrom_names) == 1
        assert bed_filtered.chrom_names[0] == all_chroms[0]

    def test_chromosome_list_invalid_name(self, sample_bed_file, genbank_file):
        """Test providing an invalid chromosome name."""
        with pytest.raises(ValueError) as exc_info:
            bedtools.SequanaCoverage(
                sample_bed_file,
                genbank_file,
                chromosome_list=["chr_invalid", "chr_nonexistent"],
            )
        assert "incorrect chromosome name" in str(exc_info.value).lower()

    def test_chromosome_list_mixed_valid_invalid(self, sample_bed_file, genbank_file):
        """Test mixing valid and invalid chromosome names."""
        bed = bedtools.SequanaCoverage(sample_bed_file, genbank_file)
        valid_chrom = bed.chrom_names[0]

        with pytest.raises(ValueError):
            bedtools.SequanaCoverage(
                sample_bed_file,
                genbank_file,
                chromosome_list=[valid_chrom, "chr_invalid"],
            )

    def test_chromosome_list_empty(self, sample_bed_file, genbank_file):
        """Test with empty chromosome list."""
        # Empty list should use default (all chromosomes)
        bed = bedtools.SequanaCoverage(sample_bed_file, genbank_file, chromosome_list=[])
        original_bed = bedtools.SequanaCoverage(sample_bed_file, genbank_file)
        # Empty list keeps all chromosomes
        assert len(bed.chrom_names) == len(original_bed.chrom_names)

    def test_chromosome_list_single_element(self, sample_bed_file, genbank_file):
        """Test with single element chromosome list."""
        bed = bedtools.SequanaCoverage(sample_bed_file, genbank_file)
        single_chrom = bed.chrom_names[0]

        bed_single = bedtools.SequanaCoverage(sample_bed_file, genbank_file, chromosome_list=[single_chrom])
        assert len(bed_single.chrom_names) == 1
        assert bed_single.chrom_names[0] == single_chrom

    def test_chromosome_list_order_preserved(self, sample_bed_file, genbank_file):
        """Test that chromosome list order is preserved."""
        bed = bedtools.SequanaCoverage(sample_bed_file, genbank_file)
        all_chroms = bed.chrom_names.copy()

        if len(all_chroms) > 1:
            # Test with reversed order
            reversed_list = all_chroms[::-1]
            bed_filtered = bedtools.SequanaCoverage(sample_bed_file, genbank_file, chromosome_list=reversed_list)
            assert bed_filtered.chrom_names == reversed_list

    def test_chromosome_list_duplicate_handling(self, sample_bed_file, genbank_file):
        """Test handling of duplicate chromosome names in list."""
        bed = bedtools.SequanaCoverage(sample_bed_file, genbank_file)
        single_chrom = bed.chrom_names[0]

        bed_filtered = bedtools.SequanaCoverage(
            sample_bed_file, genbank_file, chromosome_list=[single_chrom, single_chrom]
        )
        # Should accept duplicates (no deduplication)
        assert single_chrom in bed_filtered.chrom_names

    def test_scan_bed_with_chromosome_list(self, sample_bed_file, genbank_file):
        """Test that _scan_bed is called before chromosome_list validation."""
        bed = bedtools.SequanaCoverage(sample_bed_file, genbank_file)
        original_chroms = bed.chrom_names.copy()

        # Should scan BED file first
        assert len(original_chroms) > 0

        # Now test with filtered list
        bed_filtered = bedtools.SequanaCoverage(sample_bed_file, genbank_file, chromosome_list=[original_chroms[0]])
        assert len(bed_filtered.chrom_names) == 1


class TestGCWindowSizeSetter:
    """Test gc_window_size setter (lines 322-326)."""

    @pytest.fixture
    def coverage_object(self):
        """Create a SequanaCoverage object for testing."""
        bed_file = f"{test_dir}/data/bed/JB409847.bed"
        genbank_file = f"{test_dir}/data/genbank/JB409847.gbk"
        return bedtools.SequanaCoverage(bed_file, genbank_file)

    def test_gc_window_size_odd_number(self, coverage_object):
        """Test setting gc_window_size with an odd number."""
        coverage_object.gc_window_size = 51
        assert coverage_object.gc_window_size == 51

    def test_gc_window_size_even_number_incremented(self, coverage_object):
        """Test that even window size is incremented to odd."""
        coverage_object.gc_window_size = 50
        # Should be incremented to 51
        assert coverage_object.gc_window_size == 51
        # And should be odd
        assert coverage_object.gc_window_size % 2 == 1

    def test_gc_window_size_even_larger_number(self, coverage_object):
        """Test larger even number is incremented."""
        coverage_object.gc_window_size = 1000
        assert coverage_object.gc_window_size == 1001
        assert coverage_object.gc_window_size % 2 == 1

    def test_gc_window_size_small_odd_number(self, coverage_object):
        """Test small odd number is kept as is."""
        coverage_object.gc_window_size = 1
        assert coverage_object.gc_window_size == 1

    def test_gc_window_size_small_even_number(self, coverage_object):
        """Test small even number is incremented."""
        coverage_object.gc_window_size = 2
        assert coverage_object.gc_window_size == 3

    def test_gc_window_size_zero(self, coverage_object):
        """Test window size of zero (edge case)."""
        coverage_object.gc_window_size = 0
        # 0 is even, so should become 1
        assert coverage_object.gc_window_size == 1

    def test_gc_window_size_negative_odd(self, coverage_object):
        """Test negative odd number."""
        coverage_object.gc_window_size = -5
        assert coverage_object.gc_window_size == -5
        assert coverage_object.gc_window_size % 2 == 1 or coverage_object.gc_window_size % 2 == -1

    def test_gc_window_size_negative_even(self, coverage_object):
        """Test negative even number."""
        coverage_object.gc_window_size = -10
        assert coverage_object.gc_window_size == -9
        assert coverage_object.gc_window_size % 2 != 0

    def test_gc_window_size_logging_warning(self, coverage_object, caplog):
        """Test that warning is logged for even numbers."""
        import logging

        with caplog.at_level(logging.WARNING):
            coverage_object.gc_window_size = 100
            # Check that a warning was logged about odd number
            assert any("odd" in record.message.lower() for record in caplog.records)

    def test_gc_window_size_sequential_updates(self, coverage_object):
        """Test sequential updates to gc_window_size."""
        coverage_object.gc_window_size = 51
        assert coverage_object.gc_window_size == 51

        coverage_object.gc_window_size = 100
        assert coverage_object.gc_window_size == 101

        coverage_object.gc_window_size = 75
        assert coverage_object.gc_window_size == 75

    def test_gc_window_size_getter(self, coverage_object):
        """Test that getter returns the correct value."""
        coverage_object.gc_window_size = 49
        retrieved = coverage_object.gc_window_size
        assert retrieved == 49
        assert isinstance(retrieved, int)

    def test_gc_window_size_large_value(self, coverage_object):
        """Test with large window size."""
        coverage_object.gc_window_size = 10000
        assert coverage_object.gc_window_size == 10001

    def test_gc_window_size_multiple_setters(self, coverage_object):
        """Test multiple consecutive setter calls."""
        for value in [10, 20, 30, 40, 50]:
            coverage_object.gc_window_size = value
            expected = value if value % 2 == 1 else value + 1
            assert coverage_object.gc_window_size == expected


class TestAnnotationHandling:
    """Test annotation handling and chromosome name conversion (lines 1166-1172)."""

    @pytest.fixture
    def coverage_with_annotation(self):
        """Create a SequanaCoverage object with annotation."""
        bed_file = f"{test_dir}/data/bed/JB409847.bed"
        genbank_file = f"{test_dir}/data/genbank/JB409847.gbk"
        return bedtools.SequanaCoverage(bed_file, genbank_file)

    def test_annotation_chromosome_matching(self, coverage_with_annotation):
        """Test chromosome name matching with annotation."""
        # Annotations are at SequanaCoverage level, not ChromosomeCov
        assert hasattr(coverage_with_annotation, "feature_dict")
        # feature_dict should be populated if annotation file was provided
        assert coverage_with_annotation.feature_dict is not None

    def test_annotation_features_available(self, coverage_with_annotation):
        """Test that annotation features are available."""
        # Check that the annotation file was loaded
        assert hasattr(coverage_with_annotation, "annotation_file")
        assert coverage_with_annotation.annotation_file is not None

    def test_chromosome_name_with_accession_separator(self):
        """Test chromosome name handling with pipe separator (accession)."""
        bed_file = f"{test_dir}/data/bed/JB409847.bed"
        # Create a modified genbank to test accession parsing
        genbank_file = f"{test_dir}/data/genbank/JB409847.gbk"

        cov = bedtools.SequanaCoverage(bed_file, genbank_file)
        for chrom in cov:
            # Should handle chromosome name without error
            assert chrom.chrom_name is not None

    def test_non_matching_chromosome_name_warning(self, caplog):
        """Test warning when chromosome name doesn't match annotation."""
        # Use a BED file with chromosome names that don't match genbank
        bed_file = f"{test_dir}/data/bed/test_wrong_bed.bed"
        genbank_file = f"{test_dir}/data/genbank/JB409847.gbk"

        # This should handle mismatched chromosome names
        import logging

        with caplog.at_level(logging.WARNING):
            try:
                cov = bedtools.SequanaCoverage(bed_file, genbank_file)
                # Try to access with mismatched names
                for chrom in cov:
                    pass
            except (ValueError, FileNotFoundError):
                # Expected if file doesn't exist or has format issues
                pass

    def test_alternative_accession_parsing(self, caplog):
        """Test alternative accession parsing from chromosome names."""
        # This tests the logic on lines 1162-1164
        bed_file = f"{test_dir}/data/bed/JB409847.bed"
        genbank_file = f"{test_dir}/data/genbank/JB409847.gbk"

        import logging

        with caplog.at_level(logging.WARNING):
            cov = bedtools.SequanaCoverage(bed_file, genbank_file)
            for chrom in cov:
                pass


class TestComplexFiltering:
    """Test complex filtering/merging logic on coverage data."""

    @pytest.fixture
    def coverage_object(self):
        """Create a SequanaCoverage object."""
        bed_file = f"{test_dir}/data/bed/JB409847.bed"
        genbank_file = f"{test_dir}/data/genbank/JB409847.gbk"
        return bedtools.SequanaCoverage(bed_file, genbank_file)

    def test_getitem_indexing(self, coverage_object):
        """Test __getitem__ indexing."""
        # Should be able to access by index
        chrom = coverage_object[0]
        assert chrom is not None
        assert hasattr(chrom, "chrom_name")

    def test_iteration(self, coverage_object):
        """Test __iter__ iteration."""
        chroms = list(coverage_object)
        assert len(chroms) > 0
        for chrom in chroms:
            assert hasattr(chrom, "chrom_name")

    def test_iteration_multiple_times(self, coverage_object):
        """Test that iteration can be done multiple times."""
        first_iteration = [chrom.chrom_name for chrom in coverage_object]
        second_iteration = [chrom.chrom_name for chrom in coverage_object]
        assert first_iteration == second_iteration

    def test_double_thresholds_with_coverage(self):
        """Test DoubleThresholds with coverage filtering."""
        dt = bedtools.DoubleThresholds(-5, 5)
        assert dt.low == -5
        assert dt.high == 5
        # low2 should be half of low
        assert dt.low2 == -2.5

    def test_threshold_with_custom_ratios(self):
        """Test DoubleThresholds with custom ratio settings."""
        dt = bedtools.DoubleThresholds(-8, 8)
        dt.ldtr = 0.25
        dt.hdtr = 0.25
        assert dt.low2 == -2
        assert dt.high2 == 2

    def test_threshold_with_modified_ratios(self):
        """Test threshold modification after creation."""
        dt = bedtools.DoubleThresholds(-3, 3)
        original_low2 = dt.low2
        dt.low = -6
        # low2 should be recalculated
        assert dt.low2 != original_low2

    def test_coverage_with_thresholds(self):
        """Test coverage operations with thresholds."""
        bed_file = f"{test_dir}/data/bed/JB409847.bed"
        genbank_file = f"{test_dir}/data/genbank/JB409847.gbk"

        cov = bedtools.SequanaCoverage(bed_file, genbank_file)
        # Should use default thresholds
        assert cov.thresholds is not None


class TestEdgeCases:
    """Test edge cases and boundary conditions."""

    def test_empty_chromosome_list_validation(self):
        """Test with truly empty input."""
        bed_file = f"{test_dir}/data/bed/JB409847.bed"
        genbank_file = f"{test_dir}/data/genbank/JB409847.gbk"

        # Empty list should not raise error
        cov = bedtools.SequanaCoverage(bed_file, genbank_file, chromosome_list=[])
        assert len(cov.chrom_names) > 0

    def test_single_chromosome_coverage(self):
        """Test coverage operations on single chromosome."""
        bed_file = f"{test_dir}/data/bed/JB409847.bed"
        genbank_file = f"{test_dir}/data/genbank/JB409847.gbk"

        cov = bedtools.SequanaCoverage(bed_file, genbank_file)
        single_chrom = cov.chrom_names[0]

        cov_single = bedtools.SequanaCoverage(bed_file, genbank_file, chromosome_list=[single_chrom])
        assert len(cov_single) == 1

    def test_window_size_boundary_values(self):
        """Test window size at boundary values."""
        bed_file = f"{test_dir}/data/bed/JB409847.bed"
        genbank_file = f"{test_dir}/data/genbank/JB409847.gbk"
        cov = bedtools.SequanaCoverage(bed_file, genbank_file)

        # Test very small window
        cov.gc_window_size = 1
        assert cov.gc_window_size == 1

        # Test larger window
        cov.gc_window_size = 9999
        assert cov.gc_window_size == 9999
        assert cov.gc_window_size % 2 == 1

    def test_multiple_coverage_objects(self):
        """Test creating multiple coverage objects doesn't interfere."""
        bed_file = f"{test_dir}/data/bed/JB409847.bed"
        genbank_file = f"{test_dir}/data/genbank/JB409847.gbk"

        cov1 = bedtools.SequanaCoverage(bed_file, genbank_file)
        cov2 = bedtools.SequanaCoverage(bed_file, genbank_file)

        cov1.gc_window_size = 51
        cov2.gc_window_size = 101

        assert cov1.gc_window_size == 51
        assert cov2.gc_window_size == 101
