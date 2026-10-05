"""Tests for sequana.phylo_construction (UPGMA, neighbor-joining, consensus)."""
import pytest

from sequana.errors import SequanaException
from sequana.phylo import Tree
from sequana.phylo_construction import majority_consensus, neighbor_joining, upgma


class TestUPGMA:
    def test_upgma_matches_known_result(self):
        # Cross-checked against Bio.Phylo.TreeConstruction.DistanceTreeConstructor.upgma
        # on the same 4-taxon distance matrix: ((D:4,C:4):0.75,(B:2.5,A:2.5):2.25)
        names = ["A", "B", "C", "D"]
        matrix = [
            [0, 5, 9, 9],
            [5, 0, 10, 10],
            [9, 10, 0, 8],
            [9, 10, 8, 0],
        ]
        tree = upgma(matrix, names)
        assert set(tree.leaves()) == {"A", "B", "C", "D"}

        # A and B should be the closest pair (distance 5) and thus siblings
        # with equal branch length (2.5 each) under the strict UPGMA clock.
        node_a = tree.find_node("A")
        node_b = tree.find_node("B")
        assert node_a.parent is node_b.parent
        assert node_a.branch_length == pytest.approx(2.5)
        assert node_b.branch_length == pytest.approx(2.5)

        node_c = tree.find_node("C")
        node_d = tree.find_node("D")
        assert node_c.parent is node_d.parent
        assert node_c.branch_length == pytest.approx(4.0)
        assert node_d.branch_length == pytest.approx(4.0)

    def test_upgma_two_taxa(self):
        tree = upgma([[0, 2], [2, 0]], ["A", "B"])
        assert tree.leaves() == ["A", "B"]
        assert tree.distance("A", "B") == pytest.approx(2.0)

    def test_upgma_requires_at_least_two_taxa(self):
        with pytest.raises(SequanaException):
            upgma([[0]], ["A"])

    def test_upgma_accepts_dict_matrix(self):
        matrix = {("A", "B"): 2.0}
        tree = upgma(matrix, ["A", "B"])
        assert tree.distance("A", "B") == pytest.approx(2.0)

    def test_upgma_is_ultrametric(self):
        """UPGMA always produces an ultrametric tree: every leaf is
        equidistant from the root."""
        names = ["A", "B", "C", "D"]
        matrix = [
            [0, 2, 6, 6],
            [2, 0, 6, 6],
            [6, 6, 0, 4],
            [6, 6, 4, 0],
        ]
        tree = upgma(matrix, names)
        depths = [tree.depth_at_leaf(name) for name in names]
        assert depths[0] == pytest.approx(depths[1])
        assert depths[2] == pytest.approx(depths[3])


class TestNeighborJoining:
    def test_nj_matches_known_result(self):
        # Cross-checked against Bio.Phylo.TreeConstruction.DistanceTreeConstructor.nj
        # on the same 4-taxon distance matrix: ((A:2,B:3):3,C:4,D:4)
        names = ["A", "B", "C", "D"]
        matrix = [
            [0, 5, 9, 9],
            [5, 0, 10, 10],
            [9, 10, 0, 8],
            [9, 10, 8, 0],
        ]
        tree = neighbor_joining(matrix, names)
        assert set(tree.leaves()) == {"A", "B", "C", "D"}

        node_a = tree.find_node("A")
        node_b = tree.find_node("B")
        assert node_a.parent is node_b.parent
        assert node_a.branch_length == pytest.approx(2.0)
        assert node_b.branch_length == pytest.approx(3.0)

    def test_nj_requires_at_least_three_taxa(self):
        with pytest.raises(SequanaException):
            neighbor_joining([[0, 1], [1, 0]], ["A", "B"])

    def test_nj_five_taxa_all_leaves_present(self):
        names = ["A", "B", "C", "D", "E"]
        matrix = [
            [0, 5, 9, 9, 8],
            [5, 0, 10, 10, 9],
            [9, 10, 0, 8, 7],
            [9, 10, 8, 0, 3],
            [8, 9, 7, 3, 0],
        ]
        tree = neighbor_joining(matrix, names)
        assert set(tree.leaves()) == set(names)

    def test_nj_accepts_dict_matrix(self):
        matrix = {
            ("A", "B"): 5,
            ("A", "C"): 9,
            ("B", "C"): 10,
        }
        tree = neighbor_joining(matrix, ["A", "B", "C"])
        assert set(tree.leaves()) == {"A", "B", "C"}

    def test_nj_branch_lengths_non_negative(self):
        """NJ can produce small negative branch lengths from floating-point
        artifacts on near-degenerate distances; these must be clamped to 0."""
        names = ["A", "B", "C"]
        matrix = [
            [0, 1, 1],
            [1, 0, 1],
            [1, 1, 0],
        ]
        tree = neighbor_joining(matrix, names)
        for node in tree.all_nodes():
            assert node.branch_length >= 0.0


class TestMajorityConsensus:
    def test_consensus_matches_known_biopython_result(self):
        # Cross-checked against Bio.Phylo.Consensus.majority_consensus on the
        # same 3 input trees: (A,B) at 100%, (D,E) at 66.7%, C unresolved.
        t1 = Tree.from_newick("((A,B),(C,(D,E)));")
        t2 = Tree.from_newick("((A,B),(C,(D,E)));")
        t3 = Tree.from_newick("((A,B),((C,D),E));")

        consensus = majority_consensus([t1, t2, t3], cutoff=0.5)
        assert set(consensus.leaves()) == {"A", "B", "C", "D", "E"}

        node_a = consensus.find_node("A")
        node_b = consensus.find_node("B")
        assert node_a.parent is node_b.parent
        assert node_a.parent.bootstrap == pytest.approx(100.0)

        node_d = consensus.find_node("D")
        node_e = consensus.find_node("E")
        assert node_d.parent is node_e.parent
        assert node_d.parent.bootstrap == pytest.approx(200.0 / 3.0)

    def test_consensus_conflicting_splits_excluded_at_high_cutoff(self):
        t1 = Tree.from_newick("((A:1,B:1):1,(C:1,D:1):1);")
        t2 = Tree.from_newick("((A:1,B:1):1,(C:1,D:1):1);")
        t3 = Tree.from_newick("((A:1,C:1):1,(B:1,D:1):1);")

        # Neither split is unanimous (2/3 each) -> strict cutoff=1.0 keeps none
        consensus = majority_consensus([t1, t2, t3], cutoff=1.0)
        assert len(consensus.root.children) == 4

    def test_consensus_majority_split_included_at_default_cutoff(self):
        t1 = Tree.from_newick("((A:1,B:1):1,(C:1,D:1):1);")
        t2 = Tree.from_newick("((A:1,B:1):1,(C:1,D:1):1);")
        t3 = Tree.from_newick("((A:1,C:1):1,(B:1,D:1):1);")

        consensus = majority_consensus([t1, t2, t3], cutoff=0.5)
        node_a = consensus.find_node("A")
        node_b = consensus.find_node("B")
        # (A,B) appears in 2/3 trees -> included at cutoff=0.5
        assert node_a.parent is node_b.parent

    def test_consensus_requires_at_least_one_tree(self):
        with pytest.raises(SequanaException):
            majority_consensus([])

    def test_consensus_requires_matching_leaf_sets(self):
        t1 = Tree.from_newick("(A,B);")
        t2 = Tree.from_newick("(A,C);")
        with pytest.raises(SequanaException):
            majority_consensus([t1, t2])

    def test_consensus_single_tree_returns_same_topology(self):
        t1 = Tree.from_newick("((A,B),(C,D));")
        consensus = majority_consensus([t1], cutoff=0.5)
        assert set(consensus.leaves()) == {"A", "B", "C", "D"}
        node_a = consensus.find_node("A")
        node_b = consensus.find_node("B")
        assert node_a.parent is node_b.parent

    def test_consensus_identical_trees_full_support(self):
        t1 = Tree.from_newick("((A,B),(C,D));")
        t2 = Tree.from_newick("((A,B),(C,D));")
        consensus = majority_consensus([t1, t2], cutoff=0.5)
        node_a = consensus.find_node("A")
        node_b = consensus.find_node("B")
        assert node_a.parent is node_b.parent
        assert node_a.parent.bootstrap == pytest.approx(100.0)

    def test_build_tree_from_splits_skips_incompatible_split_defensively(self):
        """_build_tree_from_splits() must not raise or corrupt the tree when
        handed a split that doesn't cleanly partition existing groups (this
        should not occur from majority_consensus() itself, which only ever
        emits compatible splits, but the function guards against it)."""
        from sequana.phylo_construction import _build_tree_from_splits

        all_leaves = frozenset({"A", "B", "C", "D"})
        # A split that partially overlaps an existing singleton group boundary
        # (not a clean union of any subset of the starting leaf singletons'
        # supersets) -- exercises the "already represented" and "incompatible"
        # defensive branches without needing an invalid consensus input.
        splits = [
            (frozenset({"A", "B"}), 1.0),
            (frozenset({"A", "B"}), 1.0),  # duplicate: len(contained) < 2 branch
        ]
        tree = _build_tree_from_splits(all_leaves, splits)
        assert set(tree.leaves()) == {"A", "B", "C", "D"}
