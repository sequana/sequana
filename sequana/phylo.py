#
#  This file is part of Sequana software
#
#  Copyright (c) 2025 - Sequana Development Team
#
#  Distributed under the terms of the 3-clause BSD license.
#  The full license is in the LICENSE file, distributed with this software.
#
#  website: https://github.com/sequana/sequana
#  documentation: http://sequana.readthedocs.io
#
##############################################################################
"""Phylogenetic tree parsing and manipulation."""

import re
from collections import defaultdict
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Set, Tuple

import colorlog

from sequana.lazy import numpy as np

logger = colorlog.getLogger(__name__)

__all__ = ["Tree", "TreeNode"]


@dataclass
class TreeNode:
    """Single node in phylogenetic tree.

    Attributes:
        name: taxon name (leaf) or internal node label
        branch_length: distance to parent
        bootstrap: bootstrap/confidence value (0-100)
        metadata: arbitrary node annotations
    """

    name: Optional[str] = None
    branch_length: float = 0.0
    bootstrap: Optional[float] = None
    children: List["TreeNode"] = field(default_factory=list)
    parent: Optional["TreeNode"] = None
    metadata: Dict = field(default_factory=dict)

    def is_leaf(self) -> bool:
        """Return True if node has no children."""
        return len(self.children) == 0

    def is_root(self) -> bool:
        """Return True if node has no parent."""
        return self.parent is None

    def add_child(self, child: "TreeNode") -> None:
        """Add child and set its parent reference."""
        self.children.append(child)
        child.parent = self

    def __repr__(self) -> str:
        label = self.name or f"internal_{id(self)}"
        if self.bootstrap is not None:
            label += f"({self.bootstrap:.0f})"
        return label


class Tree:
    """Phylogenetic tree from Newick format.

    Parses and manipulates phylogenetic trees. Supports:
    - Standard Newick: (A:1.0,B:1.0)C:0.0;
    - Bootstrap: (A:1.0,B:1.0)95:0.0;
    - Named internal nodes: (A:1.0,B:1.0)AB:0.0;

    Examples:
        >>> tree = Tree.from_newick("(A:1.0,B:1.0)root:0.0;")
        >>> tree.leaves()
        ['A', 'B']
        >>> tree.distance("A", "B")
        2.0
    """

    def __init__(self, root):
        """Initialize tree from root node, filename, or Newick string.

        Args:
            root: TreeNode object, filename (str), or Newick format string (str)
        """
        import pathlib

        if isinstance(root, TreeNode):
            self.root = root
        elif isinstance(root, str):
            # Try to load as file first
            if pathlib.Path(root).exists():
                with open(root) as f:
                    newick_str = f.read().strip()
                parsed_root, _ = self._parse_newick(newick_str.rstrip(";"), 0)
                if parsed_root is None:
                    raise ValueError(f"Invalid Newick format in file: {root}")
                self.root = parsed_root
            else:
                # Try to parse as Newick string
                parsed_root, _ = self._parse_newick(root.rstrip(";"), 0)
                if parsed_root is None:
                    raise ValueError(f"Invalid input: not a file or valid Newick format: {root}")
                self.root = parsed_root
        else:
            raise TypeError(f"root must be TreeNode, filename, or Newick string, got {type(root)}")

        self._leaf_cache = None
        self._node_cache = None

    @classmethod
    def from_newick(cls, newick_str: str) -> "Tree":
        """Parse Newick format string.

        Format: (child1:branch1,child2:branch2)parent:branch;

        Args:
            newick_str: Newick string (with or without trailing ;)

        Returns:
            Tree object
        """
        newick_str = newick_str.strip()
        if newick_str.endswith(";"):
            newick_str = newick_str[:-1]

        root, _ = cls._parse_newick(newick_str, 0)
        if root is None:
            raise ValueError(f"Invalid Newick format: {newick_str}")
        return cls(root)

    @classmethod
    def _parse_newick(cls, s: str, pos: int) -> Tuple[Optional[TreeNode], int]:
        """Recursive Newick parser."""
        node = TreeNode()

        if pos >= len(s):
            return None, pos

        # Leaf node (name only) or internal node (recursion starts with '(')
        if s[pos] == "(":
            pos += 1
            while pos < len(s) and s[pos] != ")":
                child, pos = cls._parse_newick(s, pos)
                if child:
                    node.add_child(child)
                if pos < len(s) and s[pos] == ",":
                    pos += 1

            if pos >= len(s) or s[pos] != ")":
                return None, pos
            pos += 1

        # Parse name (if present)
        name_match = re.match(r"([^():,;]+)", s[pos:])
        if name_match:
            name = name_match.group(1)
            # Check if it's a bootstrap value (numeric) or name
            try:
                node.bootstrap = float(name)
                node.name = None
            except ValueError:
                node.name = name
            pos += len(name)

        # Parse branch length
        if pos < len(s) and s[pos] == ":":
            pos += 1
            branch_match = re.match(r"([0-9.e+-]+)", s[pos:])
            if branch_match:
                node.branch_length = float(branch_match.group(1))
                pos += len(branch_match.group(1))

        return node, pos

    def leaves(self) -> List[str]:
        """Return sorted list of leaf names."""
        if self._leaf_cache is None:
            leaves = []
            for node in self._postorder(self.root):
                if node.is_leaf() and node.name:
                    leaves.append(node.name)
            self._leaf_cache = sorted(leaves)
        return self._leaf_cache

    def leaf_count(self) -> int:
        """Return number of leaf nodes."""
        return len(self.leaves())

    def all_nodes(self) -> List[TreeNode]:
        """Return all nodes in tree (postorder traversal)."""
        if self._node_cache is None:
            self._node_cache = list(self._postorder(self.root))
        return self._node_cache

    def find_node(self, name: str) -> Optional[TreeNode]:
        """Find node by name."""
        for node in self.all_nodes():
            if node.name == name:
                return node
        return None

    def distance(self, leaf1: str, leaf2: str) -> float:
        """Euclidean distance between two leaves.

        Sum of branch lengths along path between leaves.
        """
        node1 = self.find_node(leaf1)
        node2 = self.find_node(leaf2)
        if not node1 or not node2:
            raise ValueError(f"Leaf not found: {leaf1 if not node1 else leaf2}")

        lca = self._lowest_common_ancestor(node1, node2)
        return self._distance_to_node(node1, lca) + self._distance_to_node(node2, lca)

    def prune(self, leaves_to_keep: Set[str]) -> "Tree":
        """Remove leaves not in set and return pruned tree.

        Removes single-child internal nodes (collapses lineage).
        """
        new_root = self._prune_node(self.root, leaves_to_keep)
        if new_root is None:
            raise ValueError("No leaves remain after pruning")
        return Tree(new_root)

    def subtree(self, mrca_leaves: Set[str]) -> "Tree":
        """Extract subtree containing given leaves (MRCA + descendants)."""
        lca = self._mrca(mrca_leaves)
        if lca is None:
            raise ValueError("No common ancestor found")
        # Disconnect from parent
        lca.parent = None
        return Tree(lca)

    def reroot(self, new_root_name: str) -> "Tree":
        """Reroot tree at specified node."""
        target = self.find_node(new_root_name)
        if not target:
            raise ValueError(f"Node not found: {new_root_name}")

        if target.is_root():
            return Tree(self.root)

        new_root = self._reroot_at(target)
        return Tree(new_root)

    def to_newick(self, include_branch_lengths: bool = True) -> str:
        """Serialize tree to Newick format."""
        return self._node_to_newick(self.root, include_branch_lengths) + ";"

    def stats(self) -> dict:
        """Return tree statistics."""
        nodes = self.all_nodes()
        leaves = self.leaves()
        internal_nodes = [n for n in nodes if not n.is_leaf()]
        return {
            "leaf_count": len(leaves),
            "internal_node_count": len(internal_nodes),
            "total_node_count": len(nodes),
            "height": self._tree_height(self.root),
            "leaves": leaves,
        }

    # Visualization

    def to_ascii(self, node: Optional[TreeNode] = None, prefix: str = "", is_last: bool = True) -> str:
        """Return ASCII tree representation (text-based).

        Example::

            root
            ├── A
            └── B
                ├── C
                └── D
        """
        if node is None:
            return self.to_ascii(self.root)

        lines = []
        current_prefix = "└── " if is_last else "├── "
        node_label = node.name or f"internal_{id(node)}"

        if node.bootstrap is not None:
            node_label += f"({node.bootstrap:.0f})"

        lines.append(prefix + current_prefix + node_label)

        if node.children:
            extension = "    " if is_last else "│   "
            for i, child in enumerate(node.children):
                is_last_child = i == len(node.children) - 1
                lines.append(self.to_ascii(child, prefix + extension, is_last_child))

        return "\n".join(lines)

    def plot_ascii(self) -> None:
        """Print ASCII tree to console."""
        print(self.to_ascii())

    def to_dict(self, node: Optional[TreeNode] = None) -> dict:
        """Convert tree to nested dict for JSON serialization.

        Returns:
            nested dict with keys: name, bootstrap, branch_length, children
        """
        if node is None:
            return self.to_dict(self.root)

        result = {
            "name": node.name or f"internal_{id(node)}",
            "branch_length": node.branch_length,
        }

        if node.bootstrap is not None:
            result["bootstrap"] = node.bootstrap

        if node.children:
            result["children"] = [self.to_dict(child) for child in node.children]

        return result

    def leaf_distances(self, leaf_name: str) -> Dict[str, float]:
        """Return distances from one leaf to all others.

        Args:
            leaf_name: starting leaf name

        Returns:
            dict mapping leaf names to distances
        """
        distances = {}
        for other_leaf in self.leaves():
            if other_leaf != leaf_name:
                distances[other_leaf] = self.distance(leaf_name, other_leaf)
        return distances

    def bifurcations(self) -> List[Tuple[List[str], List[str]]]:
        """Return all bifurcations as (left_leaves, right_leaves) tuples.

        Useful for cladogram analysis.
        """
        bifurcations = []

        def collect_bifurcations(node: Optional[TreeNode]):
            if node is None or node.is_leaf():
                return

            if len(node.children) == 2:
                left_leaves = self._get_all_leaves(node.children[0])
                right_leaves = self._get_all_leaves(node.children[1])
                bifurcations.append((left_leaves, right_leaves))

            for child in node.children:
                collect_bifurcations(child)

        collect_bifurcations(self.root)
        return bifurcations

    def _get_all_leaves(self, node: Optional[TreeNode]) -> List[str]:
        """Get all leaf names under a node."""
        if node is None:
            return []
        if node.is_leaf():
            return [node.name] if node.name else []
        leaves = []
        for child in node.children:
            leaves.extend(self._get_all_leaves(child))
        return leaves

    # Advanced visualization

    def plot_dendrogram(self, figsize: Tuple[int, int] = (12, 8)):
        """Plot tree as dendrogram using matplotlib.

        Args:
            figsize: figure size (width, height)
        """
        try:
            import matplotlib.pyplot as plt
            from scipy.cluster.hierarchy import dendrogram
            from scipy.spatial.distance import pdist, squareform
        except ImportError:
            logger.warning("matplotlib/scipy required for dendrogram. Install: pip install matplotlib scipy")
            return

        # Build distance matrix from tree
        leaves = self.leaves()
        n = len(leaves)
        dist_matrix = np.zeros((n, n))

        for i, leaf1 in enumerate(leaves):
            for j, leaf2 in enumerate(leaves):
                if i < j:
                    dist_matrix[i, j] = self.distance(leaf1, leaf2)
                    dist_matrix[j, i] = dist_matrix[i, j]

        # Convert to condensed distance matrix for dendrogram
        condensed = squareform(dist_matrix)

        # Perform hierarchical clustering
        from scipy.cluster.hierarchy import linkage

        linkage_matrix = linkage(condensed, method="average")

        # Plot
        fig, ax = plt.subplots(figsize=figsize)
        dendrogram(linkage_matrix, labels=leaves, ax=ax, leaf_rotation=90)
        ax.set_ylabel("Distance")
        ax.set_title(f"Tree: {len(leaves)} leaves")
        plt.tight_layout()
        return fig, ax

    def to_json(self) -> str:
        """Return tree as JSON string (for web visualization).

        Useful for D3, ETE, or other JavaScript tree viewers.
        """
        import json

        tree_dict = self.to_dict()
        return json.dumps(tree_dict, indent=2)

    def get_tree_balance(self) -> float:
        """Return tree balance metric (0-1, 1=perfectly balanced).

        Uses Colles-like index: ratio of actual vs max imbalance.
        Perfectly balanced binary tree = 1.0
        Completely ladder-like tree = 0.0
        """

        def subtree_leaf_count(node: Optional[TreeNode]) -> int:
            if node is None:
                return 0
            if node.is_leaf():
                return 1
            return sum(subtree_leaf_count(child) for child in node.children)

        def calculate_imbalance(node: Optional[TreeNode]) -> float:
            if node is None or len(node.children) <= 1:
                return 0.0

            # For each internal node with 2+ children, measure imbalance
            child_leaf_counts = [subtree_leaf_count(c) for c in node.children]
            if len(child_leaf_counts) == 2:
                # Binary node: imbalance = |left - right|
                imb = abs(child_leaf_counts[0] - child_leaf_counts[1])
            else:
                # Multi-node: use max deviation from mean
                mean_count = np.mean(child_leaf_counts)
                imb = max(abs(c - mean_count) for c in child_leaf_counts)

            # Recursively add imbalance from children
            return imb + sum(calculate_imbalance(c) for c in node.children)

        total_leaves = subtree_leaf_count(self.root)
        if total_leaves < 2:
            return 1.0

        actual_imbalance = calculate_imbalance(self.root)

        # Max imbalance for ladder-like tree with n leaves
        max_imbalance = sum(range(1, total_leaves))

        if max_imbalance == 0:
            return 1.0

        balance = 1.0 - (actual_imbalance / max_imbalance)
        return float(np.clip(balance, 0, 1))

    def get_tree_imbalance(self) -> float:
        """Return tree imbalance (1 - balance).

        High imbalance = unbalanced tree (like a ladder).
        """
        return 1.0 - self.get_tree_balance()

    def depth_at_leaf(self, leaf_name: str) -> float:
        """Return distance from root to leaf.

        Args:
            leaf_name: leaf sequence name

        Returns:
            cumulative branch length from root to leaf
        """
        leaf_node = self.find_node(leaf_name)
        if not leaf_node:
            return 0.0

        depth = 0.0
        current = leaf_node
        while current and not current.is_root():
            depth += current.branch_length
            current = current.parent

        return depth

    # Private helpers

    def _postorder(self, node: Optional[TreeNode]):
        """Postorder traversal (children before parent)."""
        if node is None:
            return
        for child in node.children:
            yield from self._postorder(child)
        yield node

    def _distance_to_node(self, node: TreeNode, target: TreeNode) -> float:
        """Sum branch lengths from node up to target."""
        distance = 0.0
        current = node
        while current != target and current is not None:
            distance += current.branch_length
            current = current.parent
        return distance

    def _lowest_common_ancestor(self, node1: TreeNode, node2: TreeNode) -> Optional[TreeNode]:
        """Find LCA of two nodes."""
        ancestors1 = []
        current = node1
        while current:
            ancestors1.append(current)
            current = current.parent

        current = node2
        while current:
            if current in ancestors1:
                return current
            current = current.parent
        return None

    def _mrca(self, leaf_names: Set[str]) -> Optional[TreeNode]:
        """Find MRCA of a set of leaves."""
        nodes = [node for name in leaf_names if (node := self.find_node(name)) is not None]
        if not nodes:
            return None
        if len(nodes) == 1:
            return nodes[0]

        lca = nodes[0]
        for node in nodes[1:]:
            lca = self._lowest_common_ancestor(lca, node)
            if lca is None:
                return None
        return lca

    def _prune_node(self, node: Optional[TreeNode], keep_leaves: Set[str]) -> Optional[TreeNode]:
        """Recursively prune node; return pruned subtree or None if all removed."""
        if node is None:
            return None

        if node.is_leaf():
            return node if node.name in keep_leaves else None

        new_children = []
        for child in node.children:
            pruned_child = self._prune_node(child, keep_leaves)
            if pruned_child:
                pruned_child.parent = node
                new_children.append(pruned_child)

        if not new_children:
            return None

        if len(new_children) == 1 and not node.is_root():
            # Collapse single-child node
            child = new_children[0]
            child.branch_length += node.branch_length
            return child

        node.children = new_children
        return node

    def _reroot_at(self, new_root: TreeNode) -> TreeNode:
        """Move root to new_root, reversing edges along path."""
        if new_root.is_root():
            return new_root

        # Collect path from new_root to old root
        path = []
        current = new_root
        while current.parent:
            path.append((current, current.parent))
            current = current.parent

        # Reverse edges: for each (child, parent) pair, detach child from parent and attach parent to child
        for i, (child, parent) in enumerate(reversed(path)):
            parent.children.remove(child)
            child.add_child(parent)
            # Keep original branch length for edge moving down the tree
            parent.branch_length = child.branch_length if i > 0 else 0.0

        new_root.branch_length = 0.0
        return new_root

    def _tree_height(self, node: Optional[TreeNode]) -> float:
        """Maximum distance from root to leaf."""
        if node is None or node.is_leaf():
            return 0.0
        return max(
            (self._tree_height(child) + child.branch_length for child in node.children),
            default=0.0,
        )

    def _node_to_newick(self, node: Optional[TreeNode], include_lengths: bool) -> str:
        """Recursive Newick serialization."""
        if node is None:
            return ""

        if node.is_leaf():
            result = node.name or ""
        else:
            children_str = ",".join(self._node_to_newick(child, include_lengths) for child in node.children)
            label = ""
            if node.bootstrap is not None:
                label = f"{node.bootstrap:.0f}"
            elif node.name:
                label = node.name
            result = f"({children_str}){label}"

        if include_lengths and node.branch_length > 0:
            result += f":{node.branch_length}"

        return result
