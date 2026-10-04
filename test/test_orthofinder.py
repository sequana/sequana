"""Tests for sequana.orthofinder."""
import pytest

from sequana.orthofinder import GeneTrees, OrthoFinder, OrthoGroups

ORTHOGROUPS_TSV = (
    "Orthogroup\tspeciesA\tspeciesB\n"
    "OG0000000\tgeneA1, geneA2\tgeneB1\n"
    "OG0000001\tgeneA3\tgeneB2, geneB3\n"
    "OG0000002\tgeneA4\t\n"
)

GENE_COUNT_TSV = (
    "Orthogroup\tspeciesA\tspeciesB\tTotal\n" "OG0000000\t2\t1\t3\n" "OG0000001\t1\t2\t3\n" "OG0000002\t1\t0\t1\n"
)

SINGLE_COPY_TXT = "OG0000002\n"

UNASSIGNED_TSV = "Orthogroup\tspeciesA\tspeciesB\n" "OG0000003\tgeneA5\t\n" "OG0000004\t\tgeneB4\n"


@pytest.fixture
def orthogroups_dir(tmpdir):
    tmpdir.join("Orthogroups.tsv").write(ORTHOGROUPS_TSV)
    tmpdir.join("Orthogroups.GeneCount.tsv").write(GENE_COUNT_TSV)
    tmpdir.join("Orthogroups_SingleCopyOrthologues.txt").write(SINGLE_COPY_TXT)
    tmpdir.join("Orthogroups_UnassignedGenes.tsv").write(UNASSIGNED_TSV)
    return str(tmpdir)


@pytest.fixture
def orthofinder_results_dir(tmpdir):
    results = tmpdir.mkdir("Orthogroups")
    results.join("Orthogroups.tsv").write(ORTHOGROUPS_TSV)
    results.join("Orthogroups.GeneCount.tsv").write(GENE_COUNT_TSV)
    return str(tmpdir)


class FakeGFF:
    """Minimal stand-in for sequana.gff3.GFF3, exposing only the .df attribute
    that add_annotation_Ld1S() relies on."""

    def __init__(self, rows):
        import pandas as pd

        self.df = pd.DataFrame(rows)


class TestOrthoGroups:
    def test_loads_and_merges_gene_counts(self, orthogroups_dir):
        o = OrthoGroups(orthogroups_dir)
        assert len(o.df) == 3
        assert "speciesA_size" in o.df.columns
        assert "speciesB_size" in o.df.columns
        assert o.df["Total_size"].tolist() == [3, 3, 1]

    def test_intermediate_frames_are_freed(self, orthogroups_dir):
        o = OrthoGroups(orthogroups_dir)
        assert not hasattr(o, "gene_counts")
        assert not hasattr(o, "ortho_groups")

    def test_single_copy_and_unassigned_loaded(self, orthogroups_dir):
        o = OrthoGroups(orthogroups_dir)
        assert len(o.single_copy_orthologs) == 1
        assert len(o.df_unassigned) == 2

    def test_summary_does_not_report_orthogroup_column_as_a_sample(self, orthogroups_dir, capsys):
        """Regression test: summary() used to iterate df_unassigned.columns
        without excluding the 'Orthogroup' column, printing a nonsensical
        'Found N unassigned genes in sample Orthogroup' line."""
        o = OrthoGroups(orthogroups_dir)
        o.summary()
        captured = capsys.readouterr()
        assert "sample Orthogroup" not in captured.out
        assert "sample speciesA" in captured.out
        assert "sample speciesB" in captured.out

    def test_summary_reports_correct_counts(self, orthogroups_dir, capsys):
        o = OrthoGroups(orthogroups_dir)
        o.summary()
        captured = capsys.readouterr()
        assert "Found 3 groups in sample speciesA" in captured.out
        assert "Found 2 groups in sample speciesB" in captured.out
        assert "Found 1 unassigned genes in sample speciesA" in captured.out
        assert "Found 1 unassigned genes in sample speciesB" in captured.out


class TestOrthoFinder:
    def test_loads_ortho_groups_and_gene_counts(self, orthofinder_results_dir):
        o = OrthoFinder(orthofinder_results_dir)
        assert len(o.ortho_groups) == 3
        assert "Total" in o.gene_counts.columns

    def test_missing_files_raise_file_not_found(self, tmpdir):
        with pytest.raises(FileNotFoundError):
            OrthoFinder(str(tmpdir))

    def test_summary_reports_counts(self, orthofinder_results_dir, capsys):
        o = OrthoFinder(orthofinder_results_dir)
        o.summary()
        captured = capsys.readouterr()
        assert "The orthogroup contains 3 groups" in captured.out
        assert "speciesA 3" in captured.out
        assert "speciesB 2" in captured.out

    def test_hist_gene_counts_runs(self, orthofinder_results_dir):
        import matplotlib

        matplotlib.use("Agg")
        o = OrthoFinder(orthofinder_results_dir)
        o.hist_gene_counts()  # should not raise

    def test_add_annotation_ld1s_matches_genes(self, orthofinder_results_dir):
        o = OrthoFinder(orthofinder_results_dir)
        o.ortho_groups = __import__("pandas").DataFrame({"Ld1S": ["geneX1, geneX2", "geneX3", ""]})
        gff = FakeGFF(
            {
                "genetic_type": ["gene", "gene", "gene"],
                "gene_id": ["geneX1", "geneX2", "geneX3"],
                "combinedAnnotation": ["annot1", "annot2", "annot3"],
                "seqid": ["chr1", "chr1", "chr2"],
                "start": [100, 200, 300],
            }
        )

        annotations, chromosomes, positions = o.add_annotation_Ld1S(gff)

        assert annotations == ["annot1", "annot3", ""]
        assert chromosomes == [["chr1"], ["chr2"], ""]
        assert positions == ["100 200", "300", ""]

    def test_add_annotation_ld1s_respects_custom_column(self, orthofinder_results_dir):
        """Regression test: the `column` parameter used to be accepted but
        silently ignored -- the method always read the hardcoded "Ld1S"
        column regardless of what was passed in."""
        import pandas as pd

        o = OrthoFinder(orthofinder_results_dir)
        o.ortho_groups = pd.DataFrame({"MySpecies": ["geneX1", "geneX3"]})
        gff = FakeGFF(
            {
                "genetic_type": ["gene", "gene"],
                "gene_id": ["geneX1", "geneX3"],
                "combinedAnnotation": ["annot1", "annot3"],
                "seqid": ["chr1", "chr2"],
                "start": [100, 300],
            }
        )

        annotations, chromosomes, positions = o.add_annotation_Ld1S(gff, column="MySpecies")
        assert annotations == ["annot1", "annot3"]

    def test_add_annotation_ld1s_handles_nan_cell(self, orthofinder_results_dir):
        """Regression test: a NaN cell (float, not str) used to be caught by
        a bare `except Exception`, masking any other unrelated bug the same
        way. Now narrowed to AttributeError specifically for this case."""
        import numpy as np
        import pandas as pd

        o = OrthoFinder(orthofinder_results_dir)
        o.ortho_groups = pd.DataFrame({"Ld1S": ["geneX1", np.nan, "geneX3"]})
        gff = FakeGFF(
            {
                "genetic_type": ["gene", "gene"],
                "gene_id": ["geneX1", "geneX3"],
                "combinedAnnotation": ["annot1", "annot3"],
                "seqid": ["chr1", "chr2"],
                "start": [100, 300],
            }
        )

        annotations, chromosomes, positions = o.add_annotation_Ld1S(gff)
        assert annotations == ["annot1", "", "annot3"]

    def test_add_annotation_ld1s_no_hit_returns_empty(self, orthofinder_results_dir):
        import pandas as pd

        o = OrthoFinder(orthofinder_results_dir)
        o.ortho_groups = pd.DataFrame({"Ld1S": ["unknown_gene"]})
        gff = FakeGFF(
            {
                "genetic_type": ["gene"],
                "gene_id": ["geneX1"],
                "combinedAnnotation": ["annot1"],
                "seqid": ["chr1"],
                "start": [100],
            }
        )

        annotations, chromosomes, positions = o.add_annotation_Ld1S(gff)
        assert annotations == [""]
        assert chromosomes == [""]
        assert positions == [""]


class TestGeneTrees:
    def test_run_computes_stats_per_tree(self, tmpdir):
        tmpdir.join("OG0000000_tree.txt").write("((A:0.1,B:0.2):0.05,(C:0.15,D:0.1):0.07);")
        tmpdir.join("OG0000001_tree.txt").write("(A:0.1);")

        gt = GeneTrees(str(tmpdir))
        df = gt._run()

        assert set(df.columns) == {"variances", "sizes", "mean_distances", "names"}
        assert len(df) == 2

        single_leaf_row = df[df["names"] == "OG0000001_tree.txt"].iloc[0]
        assert single_leaf_row["sizes"] == 1
        assert single_leaf_row["mean_distances"] == 0.0
        assert single_leaf_row["variances"] == 0.0

        multi_leaf_row = df[df["names"] == "OG0000000_tree.txt"].iloc[0]
        assert multi_leaf_row["sizes"] == 4
        assert multi_leaf_row["mean_distances"] > 0.0

    def test_run_no_trees_returns_empty_dataframe(self, tmpdir):
        gt = GeneTrees(str(tmpdir))
        df = gt._run()
        assert len(df) == 0


def test_orthofinder_results(orthofinder_results_dir):
    from sequana.orthofinder import OrthoFinderResults

    res = OrthoFinderResults(orthofinder_results_dir)
    assert res.species == ["speciesA", "speciesB"]

    shared = res.pairwise_shared()
    assert shared.loc[0, "Shared_Orthogroups"] == 2
    assert shared.loc[0, "Union_Orthogroups"] == 3

    inter = res.intersections()
    assert inter["count"].sum() == 3
    assert inter.loc[0, "count"] == 2  # OG0 and OG1 are in both

    assert list(res.categories()) == ["core", "core", "specific"]
    assert res.category_counts().loc["speciesA", "specific"] == 1
    assert res.copy_number_distribution().loc[1, "speciesA"] == 2

    assert list(res.expanded_families().index) == ["OG0000000", "OG0000001"]
    fam = res.expanded_families(with_genes=True)
    assert fam.loc["OG0000000", "Largest_in"] == "speciesA"
    assert fam.loc["OG0000000", "Genes"] == "geneA1, geneA2"
    spec = res.specific_orthogroups()
    assert list(spec["Genome"]) == ["speciesA"] and list(spec["N_genes"]) == [1]

    # optional files are absent
    assert res.statistics_overall() is None
    assert res.ortholog_relations() is None
    assert res.species_tree_file is None


def test_gene_tree_discordance(tmpdir):
    from sequana.orthofinder import gene_tree_discordance

    species = ["A", "B", "C_x", "D"]
    tmpdir.join("sp.nwk").write("((A:1,B:1):1,(C_x:1,D:1):1);")
    trees = tmpdir.mkdir("Gene_Trees")
    trees.join("OG1_tree.txt").write("((A_g1:1,B_g1:1):1,(C_x_g1:1,D_g1:1):1);")  # concordant
    trees.join("OG2_tree.txt").write("((A_g2:1,C_x_g2:1):1,(B_g2:1,D_g2:1):1);")  # discordant
    trees.join("OG3_tree.txt").write("((A_g3:1,A_g4:1):1,(B_g3:1,D_g3:1):1);")  # not single copy
    rf = gene_tree_discordance(str(trees), species, str(tmpdir.join("sp.nwk")))
    assert rf["OG1"] == 0 and rf["OG2"] == 1 and "OG3" not in rf
    # fewer than 4 genomes: uninformative
    assert gene_tree_discordance(str(trees), ["A", "B", "C_x"], str(tmpdir.join("sp.nwk"))).empty
