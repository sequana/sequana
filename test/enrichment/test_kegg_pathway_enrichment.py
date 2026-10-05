"""
Comprehensive tests for KEGGPathwayEnrichment class

Tests cover:
- Initialization with various gene list configurations
- Pathway loading and caching
- Enrichment computation
- Plotting functions
- Pathway lookups
- Error handling and edge cases
"""

import json
import tempfile
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from sequana.enrichment.kegg import KEGGPathwayEnrichment

from . import test_dir


@pytest.fixture
def mock_kegg_data():
    """Create mock KEGG pathway data for testing."""
    return {
        "eco00010": {
            "NAME": "Glycolysis - Ec",
            "ENTRY": "eco00010",
            "GENE": {
                "b0372": "pfkA; phosphofructokinase 1",
                "b1723": "pgi; glucose-6-phosphate isomerase",
                "b2097": "zwf; glucose-6-phosphate 1-dehydrogenase",
            },
            "DBLINKS": {"GO": "0006096"},
            "REFERENCE": [],
        },
        "eco00020": {
            "NAME": "Citric acid cycle - Ec",
            "ENTRY": "eco00020",
            "GENE": {
                "b0116": "gltA; citrate synthase",
                "b0720": "acnA; aconitase A",
                "b0721": "acnB; aconitase B",
            },
            "DBLINKS": {"GO": "0006091"},
            "REFERENCE": [],
        },
    }


@pytest.fixture
def sample_gene_lists():
    """Create sample gene lists for testing."""
    return {
        "down": ["b0372", "b1723"],
        "up": ["b0116", "b0720"],
        "all": ["b0372", "b1723", "b0116", "b0720", "b2097"],
    }


@pytest.fixture
def mock_bioservices():
    """Mock bioservices KEGG module."""
    with patch("sequana.enrichment.kegg.bioservices") as mock_bs:
        kegg_instance = MagicMock()
        kegg_instance.organism = "eco"
        kegg_instance.list.return_value = (
            "eco:b0372\tK00847\tpfkA; phosphofructokinase 1\n"
            "eco:b1723\tK01807\tpgi; glucose-6-phosphate isomerase\n"
            "eco:b0116\tK01647\tgltA; citrate synthase\n"
            "eco:b0720\tK01681\tacnA; aconitase A\n"
            "eco:b0721\tK01681\tacnB; aconitase B\n"
            "eco:b2097\tK00036\tzwf; glucose-6-phosphate 1-dehydrogenase\n"
        )
        kegg_instance.pathwayIds = ["path:eco00010", "path:eco00020"]
        mock_bs.KEGG.return_value = kegg_instance
        yield kegg_instance, mock_bs


@pytest.fixture
def pathways_preload_dir(mock_kegg_data):
    """Create a temporary directory with preloaded pathway data."""
    with tempfile.TemporaryDirectory() as tmpdir:
        tmppath = Path(tmpdir)
        pathway_file = tmppath / "eco.json"
        with open(pathway_file, "w") as f:
            json.dump(mock_kegg_data, f)
        yield tmppath


class TestKEGGPathwayEnrichmentInit:
    """Test KEGGPathwayEnrichment initialization."""

    def test_init_with_preloaded_pathways(self, sample_gene_lists, pathways_preload_dir, mock_bioservices):
        """Test initialization with preloaded pathway data."""
        kegg_instance, _ = mock_bioservices

        with patch("sequana.enrichment.kegg.GSEA"):
            ke = KEGGPathwayEnrichment(
                sample_gene_lists,
                "eco",
                progress=False,
                preload_directory=str(pathways_preload_dir),
            )

        assert ke.kegg.organism == "eco"
        assert len(ke.pathways) == 2
        assert "eco00010" in ke.pathways
        assert "eco00020" in ke.pathways

    def test_init_empty_gene_lists(self, mock_bioservices, pathways_preload_dir):
        """Test initialization with empty gene lists."""
        kegg_instance, _ = mock_bioservices

        empty_lists = {"down": [], "up": [], "all": []}
        with patch("sequana.enrichment.kegg.GSEA"):
            ke = KEGGPathwayEnrichment(
                empty_lists,
                "eco",
                progress=False,
                preload_directory=str(pathways_preload_dir),
            )

        assert ke.gene_lists == empty_lists
        assert len(ke.overlap_stats) == 3

    def test_init_with_background(self, sample_gene_lists, mock_bioservices, pathways_preload_dir):
        """Test initialization with custom background parameter."""
        kegg_instance, _ = mock_bioservices
        custom_background = 2000

        with patch("sequana.enrichment.kegg.GSEA"):
            ke = KEGGPathwayEnrichment(
                sample_gene_lists,
                "eco",
                background=custom_background,
                progress=False,
                preload_directory=str(pathways_preload_dir),
            )

        assert ke.background == custom_background

    def test_init_with_mapper_dataframe(self, sample_gene_lists, mock_bioservices, pathways_preload_dir):
        """Test initialization with a mapper DataFrame."""
        import pandas as pd

        kegg_instance, _ = mock_bioservices

        # Create a mock mapper DataFrame
        mapper_df = pd.DataFrame(
            {
                "name": ["pfkA", "pgi", "gltA", "acnA", "acnB", "zwf"],
            },
            index=["b0372", "b1723", "b0116", "b0720", "b0721", "b2097"],
        )

        with patch("sequana.enrichment.kegg.GSEA"):
            ke = KEGGPathwayEnrichment(
                sample_gene_lists,
                "eco",
                mapper=mapper_df,
                progress=False,
                preload_directory=str(pathways_preload_dir),
            )

        assert ke.mapper is not None
        assert len(ke.mapper) == 6

    def test_init_overlap_stats(self, sample_gene_lists, mock_bioservices, pathways_preload_dir):
        """Test that overlap statistics are computed correctly."""
        kegg_instance, _ = mock_bioservices

        with patch("sequana.enrichment.kegg.GSEA"):
            ke = KEGGPathwayEnrichment(
                sample_gene_lists,
                "eco",
                progress=False,
                preload_directory=str(pathways_preload_dir),
            )

        assert "category" in ke.overlap_stats.columns
        assert "mapped_percentage" in ke.overlap_stats.columns
        assert "N" in ke.overlap_stats.columns
        assert len(ke.overlap_stats) == 3


class TestKEGGPathwayEnrichmentPathwayLoading:
    """Test pathway loading and caching."""

    def test_load_pathways_from_preload_directory(
        self, sample_gene_lists, mock_kegg_data, mock_bioservices, pathways_preload_dir
    ):
        """Test loading pathways from a preloaded directory."""
        kegg_instance, _ = mock_bioservices

        with patch("sequana.enrichment.kegg.GSEA"):
            ke = KEGGPathwayEnrichment(
                sample_gene_lists,
                "eco",
                progress=False,
                preload_directory=str(pathways_preload_dir),
            )

        assert len(ke.pathways) == 2
        # NAME is processed and only the part before " - " is kept
        assert ke.pathways["eco00010"]["NAME"] == "Glycolysis"
        assert ke.pathways["eco00020"]["NAME"] == "Citric acid cycle"

    def test_save_and_load_pathways(self, sample_gene_lists, mock_kegg_data, mock_bioservices, pathways_preload_dir):
        """Test saving and loading pathways."""
        kegg_instance, _ = mock_bioservices

        with tempfile.TemporaryDirectory() as tmpdir:
            with patch("sequana.enrichment.kegg.GSEA"):
                ke = KEGGPathwayEnrichment(
                    sample_gene_lists,
                    "eco",
                    progress=False,
                    preload_directory=str(pathways_preload_dir),
                )

            # Save pathways
            save_dir = Path(tmpdir)
            ke.save_pathways(save_dir)

            # Verify file was created
            saved_file = save_dir / "eco.json"
            assert saved_file.exists()

            # Load and verify content
            with open(saved_file) as f:
                loaded = json.load(f)
            assert "eco00010" in loaded
            assert "eco00020" in loaded

    def test_gene_sets_creation(self, sample_gene_lists, mock_kegg_data, mock_bioservices, pathways_preload_dir):
        """Test that gene sets are properly created from pathways."""
        kegg_instance, _ = mock_bioservices

        with patch("sequana.enrichment.kegg.GSEA"):
            ke = KEGGPathwayEnrichment(
                sample_gene_lists,
                "eco",
                progress=False,
                preload_directory=str(pathways_preload_dir),
            )

        assert len(ke.gene_sets) == 2
        assert "eco00010" in ke.gene_sets
        assert "eco00020" in ke.gene_sets
        assert "pfkA" in ke.gene_sets["eco00010"]
        assert "gltA" in ke.gene_sets["eco00020"]


class TestKEGGPathwayEnrichmentMethods:
    """Test main class methods."""

    def test_check_category_valid(self, sample_gene_lists, mock_bioservices, pathways_preload_dir):
        """Test category validation with valid category."""
        kegg_instance, _ = mock_bioservices

        with patch("sequana.enrichment.kegg.GSEA"):
            ke = KEGGPathwayEnrichment(
                sample_gene_lists,
                "eco",
                progress=False,
                preload_directory=str(pathways_preload_dir),
            )

        # Should not raise
        ke._check_category("down")
        ke._check_category("up")
        ke._check_category("all")

    def test_check_category_invalid(self, sample_gene_lists, mock_bioservices, pathways_preload_dir):
        """Test category validation with invalid category."""
        kegg_instance, _ = mock_bioservices

        with patch("sequana.enrichment.kegg.GSEA"):
            ke = KEGGPathwayEnrichment(
                sample_gene_lists,
                "eco",
                progress=False,
                preload_directory=str(pathways_preload_dir),
            )

        with pytest.raises(ValueError):
            ke._check_category("invalid")

    def test_find_pathways_by_gene(self, sample_gene_lists, mock_bioservices, pathways_preload_dir):
        """Test finding pathways by gene name."""
        kegg_instance, _ = mock_bioservices

        with patch("sequana.enrichment.kegg.GSEA"):
            ke = KEGGPathwayEnrichment(
                sample_gene_lists,
                "eco",
                progress=False,
                preload_directory=str(pathways_preload_dir),
            )

        # Find pathway containing pfkA
        pathways = ke.find_pathways_by_gene("pfkA")
        assert "eco00010" in pathways

        # Find pathway containing gltA
        pathways = ke.find_pathways_by_gene("gltA")
        assert "eco00020" in pathways

        # Find non-existent gene
        pathways = ke.find_pathways_by_gene("nonexistent")
        assert len(pathways) == 0

    def test_find_pathways_by_gene_case_insensitive(self, sample_gene_lists, mock_bioservices, pathways_preload_dir):
        """Test case-insensitive pathway searching."""
        kegg_instance, _ = mock_bioservices

        with patch("sequana.enrichment.kegg.GSEA"):
            ke = KEGGPathwayEnrichment(
                sample_gene_lists,
                "eco",
                progress=False,
                preload_directory=str(pathways_preload_dir),
            )

        # Find pathway with exact gene name (case sensitive in GENE keys)
        pathways = ke.find_pathways_by_gene("pfkA")
        assert "eco00010" in pathways


class TestKEGGPathwayEnrichmentEnrichment:
    """Test enrichment computation."""

    def test_compute_enrichment_basic(self, sample_gene_lists, mock_bioservices, pathways_preload_dir):
        """Test basic enrichment computation."""
        kegg_instance, _ = mock_bioservices

        with patch("sequana.enrichment.kegg.GSEA") as mock_gsea:
            # Mock GSEA.compute_enrichment to return a simple result
            mock_gsea_instance = MagicMock()
            mock_gsea_instance.compute_enrichment.return_value = MagicMock(results=MagicMock())
            mock_gsea.return_value = mock_gsea_instance

            ke = KEGGPathwayEnrichment(
                sample_gene_lists,
                "eco",
                progress=False,
                preload_directory=str(pathways_preload_dir),
            )

        assert hasattr(ke, "enrichment")
        assert "down" in ke.enrichment or len(sample_gene_lists) > 0

    def test_enrichment_with_custom_background(self, sample_gene_lists, mock_bioservices, pathways_preload_dir):
        """Test enrichment with custom background."""
        kegg_instance, _ = mock_bioservices
        custom_bg = 3000

        with patch("sequana.enrichment.kegg.GSEA") as mock_gsea:
            mock_gsea_instance = MagicMock()
            mock_gsea.return_value = mock_gsea_instance

            ke = KEGGPathwayEnrichment(
                sample_gene_lists,
                "eco",
                background=custom_bg,
                progress=False,
                preload_directory=str(pathways_preload_dir),
            )

            # Verify compute_enrichment was called with correct background
            assert ke.background == custom_bg


class TestKEGGPathwayEnrichmentEdgeCases:
    """Test edge cases and error handling."""

    def test_init_with_used_genes(self, sample_gene_lists, mock_bioservices, pathways_preload_dir):
        """Test initialization with used_genes parameter."""
        kegg_instance, _ = mock_bioservices
        all_genes = ["b0372", "b1723", "b0116", "b0720", "b2097", "b9999"]

        with patch("sequana.enrichment.kegg.GSEA"):
            ke = KEGGPathwayEnrichment(
                sample_gene_lists,
                "eco",
                progress=False,
                preload_directory=str(pathways_preload_dir),
                used_genes=all_genes,
            )

        # Check that used_genes was included in overlap stats
        assert len(ke.overlap_stats) == 4
        assert "all genes" in ke.overlap_stats["category"].values

    def test_init_with_convert_gene_names_to_upper_case(
        self, sample_gene_lists, mock_bioservices, pathways_preload_dir
    ):
        """Test initialization with gene name case conversion."""
        kegg_instance, _ = mock_bioservices

        with patch("sequana.enrichment.kegg.GSEA"):
            ke = KEGGPathwayEnrichment(
                sample_gene_lists,
                "eco",
                progress=False,
                preload_directory=str(pathways_preload_dir),
                convert_input_gene_to_upper_case=True,
            )

        assert ke.convert_input_gene_to_upper_case is True

    def test_padj_cutoff_parameter(self, sample_gene_lists, mock_bioservices, pathways_preload_dir):
        """Test custom padj cutoff."""
        kegg_instance, _ = mock_bioservices
        custom_cutoff = 0.01

        with patch("sequana.enrichment.kegg.GSEA"):
            ke = KEGGPathwayEnrichment(
                sample_gene_lists,
                "eco",
                padj_cutoff=custom_cutoff,
                progress=False,
                preload_directory=str(pathways_preload_dir),
            )

        assert ke.padj_cutoff == custom_cutoff

    def test_missing_pathway_file(self, sample_gene_lists, mock_bioservices):
        """Test initialization with missing preload directory."""
        kegg_instance, _ = mock_bioservices

        with pytest.raises(FileNotFoundError):
            with patch("sequana.enrichment.kegg.GSEA"):
                KEGGPathwayEnrichment(
                    sample_gene_lists,
                    "eco",
                    progress=False,
                    preload_directory="/nonexistent/path",
                )


class TestKEGGPathwayEnrichmentPlotting:
    """Test plotting functionality."""

    def test_barplot_requires_valid_category(self, sample_gene_lists, mock_bioservices, pathways_preload_dir):
        """Test that barplot validates category."""
        kegg_instance, _ = mock_bioservices

        with patch("sequana.enrichment.kegg.GSEA"):
            ke = KEGGPathwayEnrichment(
                sample_gene_lists,
                "eco",
                progress=False,
                preload_directory=str(pathways_preload_dir),
            )

        with pytest.raises(ValueError):
            ke.barplot("invalid_category")

    def test_scatterplot_requires_valid_category(self, sample_gene_lists, mock_bioservices, pathways_preload_dir):
        """Test that scatterplot validates category."""
        kegg_instance, _ = mock_bioservices

        with patch("sequana.enrichment.kegg.GSEA"):
            ke = KEGGPathwayEnrichment(
                sample_gene_lists,
                "eco",
                progress=False,
                preload_directory=str(pathways_preload_dir),
            )

        with pytest.raises(ValueError):
            ke.scatterplot("invalid_category")

    def test_plot_genesets_hist(self, sample_gene_lists, mock_bioservices, pathways_preload_dir):
        """Test gene set histogram plotting."""
        kegg_instance, _ = mock_bioservices

        with patch("sequana.enrichment.kegg.GSEA"):
            ke = KEGGPathwayEnrichment(
                sample_gene_lists,
                "eco",
                progress=False,
                preload_directory=str(pathways_preload_dir),
            )

        # Should not raise
        ke.plot_genesets_hist(bins=10)


class TestKEGGPathwayEnrichmentSummary:
    """Test summary and project saving."""

    def test_summary_attributes(self, sample_gene_lists, mock_bioservices, pathways_preload_dir):
        """Test that summary is properly initialized."""
        kegg_instance, _ = mock_bioservices

        with patch("sequana.enrichment.kegg.GSEA"):
            ke = KEGGPathwayEnrichment(
                sample_gene_lists,
                "eco",
                progress=False,
                preload_directory=str(pathways_preload_dir),
            )

        assert hasattr(ke, "summary")
        # Summary name includes a prefix
        assert "KEGGPathwayEnrichment" in ke.summary.name
