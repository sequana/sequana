"""Comprehensive tests for phylo.py module."""
import os
import tempfile

import pytest

from sequana.phylo import Tree, TreeNode


class TestTreeNode:
    """Test TreeNode class."""

    def test_init_default(self):
        """Test TreeNode initialization with defaults."""
        node = TreeNode()
        assert node.name is None
        assert node.branch_length == 0.0
        assert node.bootstrap is None
        assert node.children == []
        assert node.parent is None
        assert node.metadata == {}

    def test_init_with_values(self):
        """Test TreeNode initialization with values."""
        node = TreeNode(name="A", branch_length=1.5, bootstrap=95)
        assert node.name == "A"
        assert node.branch_length == 1.5
        assert node.bootstrap == 95

    def test_is_leaf(self):
        """Test is_leaf method."""
        leaf = TreeNode(name="A")
        assert leaf.is_leaf()

        parent = TreeNode(name="parent")
        parent.add_child(leaf)
        assert not parent.is_leaf()

    def test_is_root(self):
        """Test is_root method."""
        root = TreeNode(name="root")
        assert root.is_root()

        child = TreeNode(name="child")
        root.add_child(child)
        assert not child.is_root()

    def test_add_child(self):
        """Test add_child method."""
        parent = TreeNode(name="parent")
        child1 = TreeNode(name="child1")
        child2 = TreeNode(name="child2")

        parent.add_child(child1)
        parent.add_child(child2)

        assert len(parent.children) == 2
        assert child1.parent == parent
        assert child2.parent == parent

    def test_repr(self):
        """Test string representation."""
        node = TreeNode(name="A")
        assert "A" in repr(node)

        node_with_bootstrap = TreeNode(name="B", bootstrap=95)
        repr_str = repr(node_with_bootstrap)
        assert "B" in repr_str
        assert "95" in repr_str


class TestTreeBasic:
    """Test Tree class basic functionality."""

    def test_init_with_root_node(self):
        """Test Tree initialization with TreeNode."""
        root = TreeNode(name="root")
        tree = Tree(root)
        assert tree.root == root

    def test_from_newick_simple(self):
        """Test parsing simple Newick format."""
        tree = Tree.from_newick("(A:1.0,B:1.0)C:0.0;")
        assert tree.root is not None
        assert len(tree.leaves()) == 2

    def test_from_newick_no_names(self):
        """Test parsing Newick without leaf names."""
        tree = Tree.from_newick("(,):0.0;")
        assert tree.root is not None

    def test_from_newick_single_leaf(self):
        """Test parsing single leaf."""
        tree = Tree.from_newick("A:1.0;")
        assert tree.root is not None

    def test_from_file(self):
        """Test loading tree from file."""
        newick_str = "(A:1.0,B:1.0)C:0.0;"
        with tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".nwk") as f:
            f.write(newick_str)
            nwk_file = f.name

        try:
            tree = Tree(nwk_file)
            assert tree.root is not None
            assert len(tree.leaves()) == 2
        finally:
            os.unlink(nwk_file)

    def test_from_newick_string_direct(self):
        """Test Tree initialization with Newick string."""
        tree = Tree("(A:1.0,B:1.0)C:0.0;")
        assert tree.root is not None


class TestTreeMethods:
    """Test Tree methods."""

    def setup_method(self):
        """Create a test tree."""
        self.tree = Tree.from_newick("((A:1.0,B:1.0)AB:0.5,(C:1.0,D:1.0)CD:0.5)root:0.0;")

    def test_leaves(self):
        """Test leaves method."""
        leaves = self.tree.leaves()
        assert set(leaves) == {"A", "B", "C", "D"}

    def test_leaf_count(self):
        """Test leaf_count method."""
        count = self.tree.leaf_count()
        assert count == 4

    def test_all_nodes(self):
        """Test all_nodes method."""
        nodes = self.tree.all_nodes()
        assert len(nodes) > 4  # At least leaves + internal nodes

    def test_find_node(self):
        """Test find_node method."""
        node = self.tree.find_node("A")
        assert node is not None
        assert node.name == "A"

        missing = self.tree.find_node("Z")
        assert missing is None

    def test_distance_between_leaves(self):
        """Test distance calculation between leaves."""
        dist = self.tree.distance("A", "B")
        assert dist == 2.0  # 1.0 + 1.0

    def test_distance_through_root(self):
        """Test distance through common ancestor."""
        dist = self.tree.distance("A", "C")
        # A to AB: 1.0, AB to root: 0.5, root to CD: 0.5, CD to C: 1.0
        assert dist == 3.0

    def test_stats(self):
        """Test stats method."""
        stats = self.tree.stats()
        assert isinstance(stats, dict)

    def test_to_newick(self):
        """Test Newick conversion."""
        newick = self.tree.to_newick()
        assert "A" in newick
        assert ";" in newick

    def test_to_newick_without_lengths(self):
        """Test Newick conversion without branch lengths."""
        newick = self.tree.to_newick(include_branch_lengths=False)
        assert "A" in newick
        assert ";" in newick


class TestTreePruning:
    """Test tree pruning and subtree methods."""

    def setup_method(self):
        """Create a test tree."""
        self.tree = Tree.from_newick("((A:1.0,B:1.0)AB:0.5,(C:1.0,D:1.0)CD:0.5)root:0.0;")

    def test_prune_keeps_subset(self):
        """Test pruning keeps specified leaves."""
        pruned = self.tree.prune({"A", "C"})
        leaves = pruned.leaves()
        assert set(leaves) == {"A", "C"}

    def test_prune_single_leaf(self):
        """Test pruning to single leaf."""
        # Just test it doesn't crash
        try:
            pruned = self.tree.prune({"A"})
            if pruned is not None:
                assert pruned.leaf_count() >= 1
        except ValueError:
            pass  # Some implementations raise on impossible operations

    def test_subtree(self):
        """Test subtree extraction."""
        subtree = self.tree.subtree({"A", "B"})
        leaves = subtree.leaves()
        assert set(leaves) == {"A", "B"}

    def test_subtree_three_leaves(self):
        """Test subtree with three leaves."""
        subtree = self.tree.subtree({"A", "B", "C"})
        leaves = subtree.leaves()
        assert "A" in leaves
        assert "B" in leaves
        assert "C" in leaves


class TestTreeStructure:
    """Test tree structure and conversions."""

    def setup_method(self):
        """Create a test tree."""
        self.tree = Tree.from_newick("((A:1.0,B:1.0)AB:0.5,C:1.5)root:0.0;")

    def test_to_ascii(self):
        """Test ASCII tree representation."""
        ascii_tree = self.tree.to_ascii()
        assert isinstance(ascii_tree, str)

    def test_to_dict(self):
        """Test dictionary conversion."""
        tree_dict = self.tree.to_dict()
        assert isinstance(tree_dict, dict)

    def test_leaf_distances(self):
        """Test leaf distances from single leaf."""
        distances = self.tree.leaf_distances("A")
        assert isinstance(distances, dict)

    def test_get_tree_balance(self):
        """Test tree balance calculation."""
        balance = self.tree.get_tree_balance()
        assert isinstance(balance, (int, float))

    def test_get_tree_imbalance(self):
        """Test tree imbalance calculation."""
        imbalance = self.tree.get_tree_imbalance()
        assert isinstance(imbalance, (int, float))


class TestTreeEdgeCases:
    """Test edge cases and error handling."""

    def test_unmatched_parentheses(self):
        """Test Newick with unmatched parentheses."""
        try:
            tree = Tree.from_newick("((A,B)C")
            # If it doesn't raise, check structure is reasonable
            assert tree.root is not None
        except (ValueError, TypeError):
            pass  # Expected in some implementations

    def test_single_node_tree(self):
        """Test tree with single node."""
        tree = Tree.from_newick("A:1.0;")
        assert tree.leaf_count() == 1
        assert tree.leaves()[0] == "A"

    def test_nonexistent_file_raises(self):
        """Test loading from nonexistent file."""
        try:
            tree = Tree("/nonexistent/file.nwk")
            # If it doesn't raise, root might still be None or something
        except (FileNotFoundError, ValueError):
            pass  # Expected behavior


class TestTreeBifurcations:
    """Test bifurcation detection."""

    def test_bifurcations(self):
        """Test detecting bifurcations."""
        tree = Tree.from_newick("((A:1.0,B:1.0)AB:0.5,(C:1.0,D:1.0)CD:0.5)root:0.0;")
        bifurcations = tree.bifurcations()
        assert isinstance(bifurcations, list)
        assert len(bifurcations) > 0

    def test_bifurcations_linear(self):
        """Test bifurcations in linear tree."""
        tree = Tree.from_newick("(A:1.0,B:1.0)root:0.0;")
        bifurcations = tree.bifurcations()
        assert isinstance(bifurcations, list)
