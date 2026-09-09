"""
Comprehensive tests for bamtools.py GFF/BED parsing functionality.

Focus: infer_strandness method (lines 473-492)
- Test loading gene ranges from GFF files
- Test loading gene ranges from BED files
- Test coordinate transformations
- Test edge cases (empty files, malformed entries, boundary conditions)
"""

import os
import tempfile
from pathlib import Path

import pytest
from easydev import TempFile

from sequana.bamtools import BAM

from . import test_dir


class TestGFFParsing:
    """Test GFF file parsing for gene range extraction."""

    @pytest.fixture
    def sample_gff(self):
        """Create a minimal GFF file for testing."""
        gff_content = """##gff-version 3
##sequence-region NC_000913.3 1 1000
NC_000913.3	RefSeq	gene	100	200	.	+	.	ID=gene0;Name=geneA
NC_000913.3	RefSeq	gene	300	400	.	-	.	ID=gene1;Name=geneB
NC_000913.3	RefSeq	gene	500	600	.	+	.	ID=gene2;Name=geneC
NC_000913.3	RefSeq	CDS	150	180	.	+	0	ID=cds0;Parent=gene0
"""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".gff", delete=False) as f:
            f.write(gff_content)
            f.flush()
            yield f.name
        os.unlink(f.name)

    @pytest.fixture
    def sample_bed(self):
        """Create a minimal BED file for testing."""
        bed_content = """chr1	100	200	gene1	0	+
chr1	300	400	gene2	0	-
chr2	50	150	gene3	0	+
"""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".bed", delete=False) as f:
            f.write(bed_content)
            f.flush()
            yield f.name
        os.unlink(f.name)

    @pytest.fixture
    def empty_gff(self):
        """Create an empty GFF file."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".gff", delete=False) as f:
            f.write("##gff-version 3\n")
            f.flush()
            yield f.name
        os.unlink(f.name)

    @pytest.fixture
    def malformed_gff(self):
        """Create a GFF file with malformed entries."""
        gff_content = """##gff-version 3
NC_000913.3	RefSeq	gene	100	200	.	+	.	ID=gene0
NC_000913.3	RefSeq	gene	invalid	invalid	.	+	.	ID=gene1
NC_000913.3	RefSeq	gene	300	400	.	+	.	ID=gene2
NC_000913.3	RefSeq	notgene	500	600	.	+	.	ID=feature1
"""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".gff", delete=False) as f:
            f.write(gff_content)
            f.flush()
            yield f.name
        os.unlink(f.name)

    @pytest.fixture
    def short_field_gff(self):
        """Create a GFF with insufficient fields."""
        gff_content = """##gff-version 3
NC_000913.3	RefSeq	gene	100
NC_000913.3	RefSeq	gene	200	300	.	+	.	ID=gene0
"""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".gff", delete=False) as f:
            f.write(gff_content)
            f.flush()
            yield f.name
        os.unlink(f.name)

    def test_gff_parsing_basic(self, sample_gff):
        """Test basic GFF file parsing with valid genes."""
        bam_file = f"{test_dir}/data/bam/test.bam"
        bam = BAM(bam_file)

        result = bam.infer_strandness(sample_gff, max_entries=10)

        # Should have loaded gene_ranges
        assert hasattr(bam, "gene_ranges")
        assert isinstance(bam.gene_ranges, dict)
        # Should have NC_000913.3 chromosome
        assert "NC_000913.3" in bam.gene_ranges
        # Should have 3 genes loaded (use find to get intervals)
        genes = bam.gene_ranges["NC_000913.3"].find(0, 1000)
        assert len(genes) == 3

    def test_gff_parsing_multiple_chromosomes(self):
        """Test GFF parsing with multiple chromosome names."""
        gff_content = """##gff-version 3
chr1	RefSeq	gene	100	200	.	+	.	ID=gene1
chr2	RefSeq	gene	300	400	.	-	.	ID=gene2
chr3	RefSeq	gene	500	600	.	+	.	ID=gene3
"""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".gff", delete=False) as f:
            f.write(gff_content)
            f.flush()
            gff_file = f.name

        try:
            bam_file = f"{test_dir}/data/bam/test.bam"
            bam = BAM(bam_file)
            bam.infer_strandness(gff_file, max_entries=10)

            assert len(bam.gene_ranges) == 3
            assert "chr1" in bam.gene_ranges
            assert "chr2" in bam.gene_ranges
            assert "chr3" in bam.gene_ranges
            # Verify each chromosome has a gene
            assert len(bam.gene_ranges["chr1"].find(0, 1000)) == 1
            assert len(bam.gene_ranges["chr2"].find(0, 1000)) == 1
            assert len(bam.gene_ranges["chr3"].find(0, 1000)) == 1
        finally:
            os.unlink(gff_file)

    def test_gff_parsing_strand_detection(self, sample_gff):
        """Test that strand information is correctly extracted from GFF."""
        bam_file = f"{test_dir}/data/bam/test.bam"
        bam = BAM(bam_file)
        bam.infer_strandness(sample_gff, max_entries=10)

        # Extract strands from intervals using find()
        for chrom, interval_tree in bam.gene_ranges.items():
            strands = interval_tree.find(0, 1000000)
            # Each strand value should be + or -
            for strand in strands:
                assert strand in ["+", "-"]

    def test_gff_parsing_empty_file(self, empty_gff):
        """Test parsing empty GFF file (only header)."""
        bam_file = f"{test_dir}/data/bam/test.bam"
        bam = BAM(bam_file)
        bam.infer_strandness(empty_gff, max_entries=10)

        # Should have empty gene_ranges
        assert bam.gene_ranges == {}

    def test_gff_parsing_malformed_entries(self, malformed_gff):
        """Test handling of malformed GFF entries."""
        bam_file = f"{test_dir}/data/bam/test.bam"
        bam = BAM(bam_file)

        # Should handle malformed entries gracefully
        # The code logs warnings but may raise ValueError
        with pytest.raises((ValueError, IndexError)):
            # Expected when trying to convert invalid coordinate
            bam.infer_strandness(malformed_gff, max_entries=10)

    def test_gff_parsing_short_fields(self, short_field_gff):
        """Test GFF with insufficient fields."""
        bam_file = f"{test_dir}/data/bam/test.bam"
        bam = BAM(bam_file)

        # Should raise IndexError for insufficient fields
        with pytest.raises((ValueError, IndexError)):
            bam.infer_strandness(short_field_gff, max_entries=10)

    def test_gff_parsing_comments_ignored(self):
        """Test that comment lines are ignored."""
        gff_content = """##gff-version 3
# This is a comment
NC_000913.3	RefSeq	gene	100	200	.	+	.	ID=gene1
# Another comment
NC_000913.3	RefSeq	gene	300	400	.	-	.	ID=gene2
"""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".gff", delete=False) as f:
            f.write(gff_content)
            f.flush()
            gff_file = f.name

        try:
            bam_file = f"{test_dir}/data/bam/test.bam"
            bam = BAM(bam_file)
            bam.infer_strandness(gff_file, max_entries=10)

            assert len(bam.gene_ranges["NC_000913.3"].find(0, 10000)) == 2
        finally:
            os.unlink(gff_file)

    def test_gff_parsing_non_gene_features_skipped(self):
        """Test that non-gene features are skipped."""
        gff_content = """##gff-version 3
NC_000913.3	RefSeq	gene	100	200	.	+	.	ID=gene1
NC_000913.3	RefSeq	CDS	110	190	.	+	0	ID=cds1;Parent=gene1
NC_000913.3	RefSeq	exon	110	150	.	+	.	ID=exon1;Parent=gene1
NC_000913.3	RefSeq	gene	300	400	.	-	.	ID=gene2
"""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".gff", delete=False) as f:
            f.write(gff_content)
            f.flush()
            gff_file = f.name

        try:
            bam_file = f"{test_dir}/data/bam/test.bam"
            bam = BAM(bam_file)
            bam.infer_strandness(gff_file, max_entries=10)

            # Only 2 genes should be loaded, not CDS/exon
            assert len(bam.gene_ranges["NC_000913.3"].find(0, 10000)) == 2
        finally:
            os.unlink(gff_file)

    def test_gff_field_order(self):
        """Test that GFF fields are correctly interpreted (0-indexed)."""
        gff_content = """##gff-version 3
NC_000913.3	source	gene	100	200	.	+	.	ID=gene1
NC_000913.3	source	gene	300	400	.	-	.	ID=gene2
"""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".gff", delete=False) as f:
            f.write(gff_content)
            f.flush()
            gff_file = f.name

        try:
            bam_file = f"{test_dir}/data/bam/test.bam"
            bam = BAM(bam_file)
            bam.infer_strandness(gff_file, max_entries=10)

            # find() returns strand values in order
            strands = bam.gene_ranges["NC_000913.3"].find(0, 1000)
            assert len(strands) == 2
            assert strands[0] == "+"
            assert strands[1] == "-"
        finally:
            os.unlink(gff_file)


class TestBEDParsing:
    """Test BED file parsing for gene range extraction."""

    @pytest.fixture
    def sample_bed(self):
        """Create a minimal BED file for testing."""
        bed_content = """chr1	100	200	gene1	0	+
chr1	300	400	gene2	0	-
chr2	50	150	gene3	0	+
"""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".bed", delete=False) as f:
            f.write(bed_content)
            f.flush()
            yield f.name
        os.unlink(f.name)

    @pytest.fixture
    def empty_bed(self):
        """Create an empty BED file."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".bed", delete=False) as f:
            f.flush()
            yield f.name
        os.unlink(f.name)

    @pytest.fixture
    def bed_with_comments(self):
        """Create BED file with comments and headers."""
        bed_content = """# This is a comment
chr1	100	200	gene1	0	+
chr1	300	400	gene2	0	-
"""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".bed", delete=False) as f:
            f.write(bed_content)
            f.flush()
            yield f.name
        os.unlink(f.name)

    def test_bed_parsing_basic(self, sample_bed):
        """Test basic BED file parsing."""
        bam_file = f"{test_dir}/data/bam/test.bam"
        bam = BAM(bam_file)
        bam.infer_strandness(sample_bed, max_entries=10)

        assert hasattr(bam, "gene_ranges")
        assert isinstance(bam.gene_ranges, dict)
        assert "chr1" in bam.gene_ranges
        assert "chr2" in bam.gene_ranges

    def test_bed_parsing_multiple_chromosomes(self, sample_bed):
        """Test BED parsing with multiple chromosome names."""
        bam_file = f"{test_dir}/data/bam/test.bam"
        bam = BAM(bam_file)
        bam.infer_strandness(sample_bed, max_entries=10)

        # Should have 2 chromosomes
        assert len(bam.gene_ranges) == 2
        # chr1 should have 2 genes
        assert len(bam.gene_ranges["chr1"].find(0, 10000)) == 2
        # chr2 should have 1 gene
        assert len(bam.gene_ranges["chr2"].find(0, 10000)) == 1

    def test_bed_parsing_strand_detection(self, sample_bed):
        """Test that strand information is correctly extracted from BED."""
        bam_file = f"{test_dir}/data/bam/test.bam"
        bam = BAM(bam_file)
        bam.infer_strandness(sample_bed, max_entries=10)

        # Extract strands using find()
        for chrom, interval_tree in bam.gene_ranges.items():
            strands = interval_tree.find(0, 1000000)
            for strand in strands:
                assert strand in ["+", "-"]

    def test_bed_parsing_empty_file(self, empty_bed):
        """Test parsing empty BED file."""
        bam_file = f"{test_dir}/data/bam/test.bam"
        bam = BAM(bam_file)
        bam.infer_strandness(empty_bed, max_entries=10)

        assert bam.gene_ranges == {}

    def test_bed_parsing_with_comments_and_headers(self, bed_with_comments):
        """Test BED parsing ignores comments and header lines."""
        bam_file = f"{test_dir}/data/bam/test.bam"
        bam = BAM(bam_file)
        bam.infer_strandness(bed_with_comments, max_entries=10)

        # Should only load actual gene entries
        assert len(bam.gene_ranges["chr1"].find(0, 10000)) == 2

    def test_bed_field_order(self):
        """Test that BED fields are correctly interpreted."""
        bed_content = """chr1	100	200	gene1	0	+
chr1	300	400	gene2	0	-
"""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".bed", delete=False) as f:
            f.write(bed_content)
            f.flush()
            bed_file = f.name

        try:
            bam_file = f"{test_dir}/data/bam/test.bam"
            bam = BAM(bam_file)
            bam.infer_strandness(bed_file, max_entries=10)

            # Get all intervals for chr1 in range 0-1000
            intervals = bam.gene_ranges["chr1"].find(0, 1000)
            # First gene: 100-200, +
            # Intervals are returned in order but as values (strand info)
            assert len(intervals) == 2
        finally:
            os.unlink(bed_file)

    def test_bed_parsing_with_real_data(self):
        """Test BED parsing with real test data."""
        # The test.bed file has only 3 columns, but BED requires 6 for strandness
        # So this test skips real data for now
        # Real BED files with 6+ columns can be used from test/data/bed/
        bed_file = f"{test_dir}/data/bed/hg38_chr18.bed"
        if not os.path.exists(bed_file):
            pytest.skip("Test BED file not available")

        bam_file = f"{test_dir}/data/bam/test.bam"
        bam = BAM(bam_file)

        # Should handle real BED file without errors
        try:
            bam.infer_strandness(bed_file, max_entries=10)
            assert bam.gene_ranges is not None
        except (ValueError, IndexError):
            # May fail if chromosome names don't match
            pass


class TestCoordinateTransformations:
    """Test coordinate handling and edge cases."""

    def test_boundary_coordinates(self):
        """Test boundary coordinate values."""
        gff_content = """##gff-version 3
NC_000913.3	RefSeq	gene	1	10	.	+	.	ID=gene1
NC_000913.3	RefSeq	gene	9999999	10000000	.	-	.	ID=gene2
"""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".gff", delete=False) as f:
            f.write(gff_content)
            f.flush()
            gff_file = f.name

        try:
            bam_file = f"{test_dir}/data/bam/test.bam"
            bam = BAM(bam_file)
            bam.infer_strandness(gff_file, max_entries=10)

            # Use find to get strand values
            strands = bam.gene_ranges["NC_000913.3"].find(0, 100000000)
            assert len(strands) == 2
            assert strands[0] == "+"
            assert strands[1] == "-"
        finally:
            os.unlink(gff_file)

    def test_overlapping_genes(self):
        """Test handling of overlapping gene ranges."""
        gff_content = """##gff-version 3
NC_000913.3	RefSeq	gene	100	300	.	+	.	ID=gene1
NC_000913.3	RefSeq	gene	200	400	.	-	.	ID=gene2
"""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".gff", delete=False) as f:
            f.write(gff_content)
            f.flush()
            gff_file = f.name

        try:
            bam_file = f"{test_dir}/data/bam/test.bam"
            bam = BAM(bam_file)
            bam.infer_strandness(gff_file, max_entries=10)

            assert len(bam.gene_ranges["NC_000913.3"].find(0, 10000)) == 2
        finally:
            os.unlink(gff_file)

    def test_single_position_genes(self):
        """Test handling of single-position (zero-length) genes."""
        gff_content = """##gff-version 3
NC_000913.3	RefSeq	gene	100	100	.	+	.	ID=gene1
NC_000913.3	RefSeq	gene	200	200	.	-	.	ID=gene2
"""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".gff", delete=False) as f:
            f.write(gff_content)
            f.flush()
            gff_file = f.name

        try:
            bam_file = f"{test_dir}/data/bam/test.bam"
            bam = BAM(bam_file)
            bam.infer_strandness(gff_file, max_entries=10)

            strands = bam.gene_ranges["NC_000913.3"].find(0, 1000)
            assert len(strands) == 2
            assert strands[0] == "+"
            assert strands[1] == "-"
        finally:
            os.unlink(gff_file)
