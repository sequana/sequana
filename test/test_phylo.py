"""Tests for phylo module."""

import tempfile
from pathlib import Path

import pytest

from sequana.phylo import Tree, TreeNode


class TestTreeNode:
    """Test TreeNode dataclass."""

    def test_is_leaf(self):
        """Leaf has no children."""
        node = TreeNode(name="A")
        assert node.is_leaf()

        child = TreeNode(name="B")
        node.add_child(child)
        assert not node.is_leaf()

    def test_is_root(self):
        """Root has no parent."""
        root = TreeNode(name="root")
        assert root.is_root()

        child = TreeNode(name="child")
        root.add_child(child)
        assert not child.is_root()
        assert root.is_root()

    def test_add_child(self):
        """Add child sets parent reference."""
        parent = TreeNode(name="parent")
        child = TreeNode(name="child")
        parent.add_child(child)

        assert child in parent.children
        assert child.parent is parent

    def test_repr(self):
        """String repr includes name and bootstrap."""
        node = TreeNode(name="A")
        assert "A" in repr(node)

        node_with_boot = TreeNode(name="B", bootstrap=95.0)
        assert "B" in repr(node_with_boot)
        assert "95" in repr(node_with_boot)


class TestTreeInit:
    """Test Tree initialization with different input types."""

    def test_init_with_treenode(self):
        """Initialize with TreeNode object (original API)."""
        root = TreeNode(name="root")
        child1 = TreeNode(name="A", branch_length=1.0)
        child2 = TreeNode(name="B", branch_length=1.0)
        root.add_child(child1)
        root.add_child(child2)

        t = Tree(root)
        assert t.root is root
        assert t.leaves() == ["A", "B"]

    def test_init_with_newick_string(self):
        """Initialize with Newick format string."""
        newick = "(A:1.0,B:1.0)root:0.0;"
        t = Tree(newick)

        assert t.root.name == "root"
        assert t.leaves() == ["A", "B"]

    def test_init_with_newick_string_no_semicolon(self):
        """Newick without trailing semicolon works."""
        newick = "(A:1.0,B:1.0)root:0.0"
        t = Tree(newick)

        assert t.root.name == "root"
        assert t.leaves() == ["A", "B"]

    def test_init_with_file(self, tmp_path):
        """Initialize by reading from file."""
        newick = "(A:1.0,B:1.0)root:0.0;"
        tree_file = tmp_path / "test.tree"
        tree_file.write_text(newick)

        t = Tree(str(tree_file))
        assert t.root.name == "root"
        assert t.leaves() == ["A", "B"]

    def test_init_with_invalid_input(self):
        """Invalid input raises TypeError."""
        with pytest.raises(TypeError):
            Tree(123)

    def test_init_with_invalid_newick(self):
        """Invalid Newick parses as single leaf (lenient parser)."""
        # Parser is lenient and parses text as leaf name
        t = Tree("not a valid newick string at all !!!")
        assert t.root is not None
        # It parses the text as a leaf name
        assert len(t.leaves()) > 0

    def test_init_with_nonexistent_file(self):
        """Nonexistent file path parsed as Newick (lenient)."""
        # Nonexistent file path is parsed as Newick string (lenient parser)
        t = Tree("/nonexistent/path/to/file.tree")
        # Parser treats it as a leaf name
        assert t.root is not None


class TestTreeFromNewick:
    """Test Tree.from_newick classmethod."""

    def test_from_newick_simple(self):
        """Parse simple two-leaf tree."""
        newick = "(A:1.0,B:1.0)root:0.0;"
        t = Tree.from_newick(newick)

        assert t.root.name == "root"
        assert t.leaves() == ["A", "B"]

    def test_from_newick_bootstrap(self):
        """Parse tree with bootstrap values."""
        newick = "(A:1.0,B:1.0)95:0.0;"
        t = Tree.from_newick(newick)

        assert t.root.bootstrap == 95.0
        assert t.leaves() == ["A", "B"]

    def test_from_newick_complex(self):
        """Parse tree with multiple levels."""
        newick = "((A:1.0,B:1.0)90:0.5,(C:0.8,D:0.8)85:0.5)root:0.0;"
        t = Tree.from_newick(newick)

        assert t.root.name == "root"
        leaves = t.leaves()
        assert set(leaves) == {"A", "B", "C", "D"}


class TestTreeMethods:
    """Test Tree methods."""

    @pytest.fixture
    def simple_tree(self):
        """(A:1.0,B:1.0)root:0.0;"""
        return Tree.from_newick("(A:1.0,B:1.0)root:0.0;")

    @pytest.fixture
    def complex_tree(self):
        """((A:1.0,B:1.0)90:0.5,(C:0.8,D:0.8)85:0.5)root:0.0;"""
        return Tree.from_newick("((A:1.0,B:1.0)90:0.5,(C:0.8,D:0.8)85:0.5)root:0.0;")

    def test_leaves(self, simple_tree):
        """Return list of leaf names."""
        assert simple_tree.leaves() == ["A", "B"]

    def test_distance(self, simple_tree):
        """Calculate distance between two leaves."""
        dist = simple_tree.distance("A", "B")
        assert dist == 2.0  # 1.0 + 1.0

    def test_distance_same_leaf(self, simple_tree):
        """Distance from leaf to itself is 0."""
        dist = simple_tree.distance("A", "A")
        assert dist == 0.0

    def test_all_nodes(self, complex_tree):
        """Return all nodes in tree."""
        nodes = complex_tree.all_nodes()
        names = [n.name for n in nodes if n.name]
        assert "A" in names
        assert "B" in names
        assert "C" in names
        assert "D" in names
        assert "root" in names

    def test_to_ascii(self, simple_tree):
        """ASCII representation."""
        ascii_repr = simple_tree.to_ascii()
        assert isinstance(ascii_repr, str)
        assert len(ascii_repr) > 0

    def test_to_newick(self, simple_tree):
        """Convert back to Newick format."""
        newick = simple_tree.to_newick()
        assert isinstance(newick, str)
        assert "A" in newick
        assert "B" in newick

    def test_to_dict(self, simple_tree):
        """Convert to dictionary structure."""
        d = simple_tree.to_dict()
        assert isinstance(d, dict)
        assert "name" in d
        assert "children" in d

    def test_to_json(self, simple_tree):
        """Convert to JSON string."""
        j = simple_tree.to_json()
        assert isinstance(j, str)
        assert "A" in j
        assert "B" in j

    def test_stats(self, simple_tree):
        """Get tree statistics."""
        stats = simple_tree.stats()
        assert isinstance(stats, dict)
        assert "leaf_count" in stats
        assert stats["leaf_count"] == 2


class TestPlotDendrogram:
    """Test plot_dendrogram method."""

    def test_plot_dendrogram_returns_fig_ax(self):
        """plot_dendrogram returns matplotlib figure and axes."""
        t = Tree.from_newick("(A:1.0,B:1.0)root:0.0;")
        fig, ax = t.plot_dendrogram()

        assert fig is not None
        assert ax is not None
        # Close to avoid display warnings
        import matplotlib.pyplot as plt

        plt.close(fig)

    def test_plot_dendrogram_figsize(self):
        """plot_dendrogram accepts figsize parameter."""
        t = Tree.from_newick("(A:1.0,B:1.0)root:0.0;")
        fig, ax = t.plot_dendrogram(figsize=(10, 6))

        assert fig.get_figwidth() == 10
        assert fig.get_figheight() == 6
        import matplotlib.pyplot as plt

        plt.close(fig)

    def test_plot_dendrogram_with_complex_tree(self):
        """plot_dendrogram works with complex tree."""
        t = Tree.from_newick("((A:1.0,B:1.0)90:0.5,(C:0.8,D:0.8)85:0.5)root:0.0;")
        fig, ax = t.plot_dendrogram()

        assert fig is not None
        import matplotlib.pyplot as plt

        plt.close(fig)


class TestRealTreeFile:
    """Test with real tree file."""

    def test_init_with_real_file(self):
        """Initialize from actual tree file if available."""
        tree_file = Path(__file__).parent / "data" / "phylo_7511scos.treefile"
        if not tree_file.exists():
            pytest.skip("Test tree file not found")

        t = Tree(str(tree_file))
        assert t.root is not None
        assert len(t.leaves()) > 0

    def test_plot_dendrogram_real_file(self):
        """Plot dendrogram from real tree file."""
        tree_file = Path(__file__).parent / "data" / "phylo_7511scos.treefile"
        if not tree_file.exists():
            pytest.skip("Test tree file not found")

        t = Tree(str(tree_file))
        fig, ax = t.plot_dendrogram()
        assert fig is not None

        import matplotlib.pyplot as plt

        plt.close(fig)
