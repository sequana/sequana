#
#  This file is part of Sequana software
#
#  Copyright (c) 2026 - Sequana Development Team
#
#  Distributed under the terms of the 3-clause BSD license.
#  The full license is in the LICENSE file, distributed with this software.
#
#  website: https://github.com/sequana/sequana
#  documentation: http://sequana.readthedocs.io
#
##############################################################################
"""Phylogenetic tree construction (UPGMA, neighbor-joining) and consensus trees.

Complements :mod:`sequana.phylo` (which parses/manipulates an already-built
``Tree``) with the missing piece: building a tree *from* a distance matrix,
and summarizing a set of bootstrap-replicate trees into a majority-rule
consensus tree. This is the Biopython ``Bio.Phylo.TreeConstruction`` /
``Bio.Phylo.Consensus`` equivalent.

Example, UPGMA from a 4-taxon distance matrix::

    from sequana.phylo_construction import upgma

    names = ["A", "B", "C", "D"]
    matrix = [
        [0, 5, 9, 9],
        [5, 0, 10, 10],
        [9, 10, 0, 8],
        [9, 10, 8, 0],
    ]
    tree = upgma(matrix, names)
    print(tree.to_newick())

Example, neighbor-joining and majority-rule consensus of bootstrap trees::

    from sequana.phylo_construction import neighbor_joining, majority_consensus

    tree = neighbor_joining(matrix, names)
    consensus = majority_consensus([tree1, tree2, tree3], cutoff=0.5)
"""
from collections import Counter
from typing import Dict, List, Optional, Sequence, Set, Tuple, Union

import colorlog

from sequana.errors import SequanaException
from sequana.phylo import Tree, TreeNode

logger = colorlog.getLogger(__name__)

__all__ = ["upgma", "neighbor_joining", "majority_consensus"]

DistanceMatrix = Union[Sequence[Sequence[float]], Dict[Tuple[str, str], float]]


def _matrix_to_dict(matrix: DistanceMatrix, names: List[str]) -> Dict[Tuple[int, int], float]:
    """Normalize a distance matrix (2D array-like or name-pair dict) to an
    ``{(i, j): distance}`` dict indexed by position in ``names``."""
    n = len(names)
    dist: Dict[Tuple[int, int], float] = {}

    if isinstance(matrix, dict):
        name_to_idx = {name: idx for idx, name in enumerate(names)}
        for (a, b), d in matrix.items():
            i, j = name_to_idx[a], name_to_idx[b]
            dist[(i, j)] = float(d)
            dist[(j, i)] = float(d)
    else:
        for i in range(n):
            for j in range(n):
                if i != j:
                    dist[(i, j)] = float(matrix[i][j])

    return dist


def upgma(matrix: DistanceMatrix, names: List[str]) -> Tree:
    """Build a rooted tree by UPGMA (Unweighted Pair Group Method with
    Arithmetic mean) clustering of a pairwise distance matrix.

    UPGMA assumes a constant molecular clock (equal evolutionary rate across
    lineages) -- it always produces an ultrametric tree (all leaves
    equidistant from the root). For divergent-rate data, prefer
    :func:`neighbor_joining`.

    Args:
        matrix: NxN distance matrix (list of lists / 2D array-like, symmetric,
            zero diagonal) or a dict mapping (name_a, name_b) -> distance.
        names: taxon names, in the same order as the matrix rows/columns.

    Returns:
        Tree: rooted tree with branch lengths set from cluster merge heights.

    Raises:
        SequanaException: if fewer than 2 names are given.

    Example::

        from sequana.phylo_construction import upgma
        tree = upgma([[0, 2], [2, 0]], ["A", "B"])
        tree.leaves()
        ['A', 'B']
    """
    n = len(names)
    if n < 2:
        raise SequanaException("upgma() requires at least 2 taxa")

    D = _matrix_to_dict(matrix, names)

    # node id -> (TreeNode, height, size)
    height: Dict[int, float] = {i: 0.0 for i in range(n)}
    size: Dict[int, int] = {i: 1 for i in range(n)}
    node: Dict[int, TreeNode] = {i: TreeNode(name=names[i]) for i in range(n)}

    active: Set[int] = set(range(n))
    next_id = n

    while len(active) > 1:
        active_list = sorted(active)
        best_pair: Optional[Tuple[int, int]] = None
        best_d: Optional[float] = None
        for a in range(len(active_list)):
            for b in range(a + 1, len(active_list)):
                i, j = active_list[a], active_list[b]
                d = D[(i, j)]
                if best_d is None or d < best_d:
                    best_d = d
                    best_pair = (i, j)

        i, j = best_pair  # type: ignore[misc]
        new_height = best_d / 2.0  # type: ignore[operator]

        merged = TreeNode()
        child_i = node[i]
        child_j = node[j]
        child_i.branch_length = new_height - height[i]
        child_j.branch_length = new_height - height[j]
        merged.add_child(child_i)
        merged.add_child(child_j)

        new_size = size[i] + size[j]
        for k in active:
            if k in (i, j):
                continue
            new_d = (size[i] * D[(i, k)] + size[j] * D[(j, k)]) / new_size
            D[(next_id, k)] = new_d
            D[(k, next_id)] = new_d

        node[next_id] = merged
        height[next_id] = new_height
        size[next_id] = new_size

        active.discard(i)
        active.discard(j)
        active.add(next_id)
        next_id += 1

    root_id = next(iter(active))
    return Tree(node[root_id])


def neighbor_joining(matrix: DistanceMatrix, names: List[str]) -> Tree:
    """Build an unrooted tree by neighbor-joining (Saitou & Nei, 1987).

    Unlike UPGMA, neighbor-joining does not assume a constant rate of
    evolution and produces the correct topology under more relaxed
    conditions (additive distances), at the cost of typically returning an
    unrooted tree topology (returned here rooted arbitrarily at the last
    internal node, matching common convention).

    Args:
        matrix: NxN distance matrix (list of lists / 2D array-like, symmetric,
            zero diagonal) or a dict mapping (name_a, name_b) -> distance.
        names: taxon names, in the same order as the matrix rows/columns.

    Returns:
        Tree: tree with branch lengths from the neighbor-joining algorithm.

    Raises:
        SequanaException: if fewer than 3 names are given (NJ needs at least
            3 taxa to have a meaningful internal topology).

    Example::

        from sequana.phylo_construction import neighbor_joining
        names = ["A", "B", "C", "D"]
        matrix = [[0, 5, 9, 9], [5, 0, 10, 10], [9, 10, 0, 8], [9, 10, 8, 0]]
        tree = neighbor_joining(matrix, names)
        tree.leaves()
        ['A', 'B', 'C', 'D']
    """
    n = len(names)
    if n < 3:
        raise SequanaException("neighbor_joining() requires at least 3 taxa")

    D = _matrix_to_dict(matrix, names)
    node: Dict[int, TreeNode] = {i: TreeNode(name=names[i]) for i in range(n)}
    active: List[int] = list(range(n))
    next_id = n

    while len(active) > 2:
        r = len(active)
        total_dist = {i: sum(D[(i, k)] for k in active if k != i) for i in active}

        best_q: Optional[float] = None
        best_pair: Optional[Tuple[int, int]] = None
        for a in range(len(active)):
            for b in range(a + 1, len(active)):
                i, j = active[a], active[b]
                q = (r - 2) * D[(i, j)] - total_dist[i] - total_dist[j]
                if best_q is None or q < best_q:
                    best_q = q
                    best_pair = (i, j)

        i, j = best_pair  # type: ignore[misc]
        dij = D[(i, j)]
        li = 0.5 * dij + (total_dist[i] - total_dist[j]) / (2 * (r - 2))
        lj = dij - li
        # Guard against tiny negative branch lengths from floating-point
        # noise or near-degenerate distances; NJ branch lengths should be
        # non-negative for a biologically meaningful tree.
        li = max(li, 0.0)
        lj = max(lj, 0.0)

        merged = TreeNode()
        child_i = node[i]
        child_j = node[j]
        child_i.branch_length = li
        child_j.branch_length = lj
        merged.add_child(child_i)
        merged.add_child(child_j)
        node[next_id] = merged

        for k in active:
            if k in (i, j):
                continue
            new_d = 0.5 * (D[(i, k)] + D[(j, k)] - dij)
            D[(next_id, k)] = new_d
            D[(k, next_id)] = new_d

        active = [x for x in active if x not in (i, j)] + [next_id]
        next_id += 1

    # Two nodes remain: join them directly with a single branch, rooting
    # (arbitrarily, as NJ trees are inherently unrooted) at the second node.
    i, j = active
    dij = max(D[(i, j)], 0.0)
    root = node[j]
    child_i = node[i]
    child_i.branch_length = dij
    root.add_child(child_i)

    return Tree(root)


def _get_bipartitions(tree: Tree, all_leaves: frozenset) -> Set[frozenset]:
    """Return the set of non-trivial bipartitions present in a tree.

    Each internal node's descendant leaf set is compared against the *whole*
    tree's leaf set (not just its sibling's leaves, unlike
    ``Tree.bifurcations()`` which only reports local child splits) -- this is
    what makes nested splits (e.g. a (D,E) clade inside a larger tree)
    comparable across different trees during consensus building.
    """
    splits: Set[frozenset] = set()
    n_leaves = len(all_leaves)

    def collect(node: Optional[TreeNode]):
        if node is None or node.is_leaf():
            return
        descendant_leaves = frozenset(tree._get_all_leaves(node))
        if 0 < len(descendant_leaves) < n_leaves:
            complement = all_leaves - descendant_leaves
            canonical = descendant_leaves if len(descendant_leaves) <= len(complement) else complement
            splits.add(canonical)
        for child in node.children:
            collect(child)

    collect(tree.root)
    return splits


def majority_consensus(trees: List[Tree], cutoff: float = 0.5) -> Tree:
    """Build a majority-rule consensus tree from a list of trees on the same taxa.

    Typical use: summarizing bootstrap replicate trees into a single
    consensus topology, with each retained internal split annotated by the
    fraction of input trees that contain it (stored as ``bootstrap`` on the
    corresponding internal node, 0-100 scale to match the existing
    ``TreeNode.bootstrap`` convention used when parsing Newick bootstrap
    values).

    Args:
        trees: list of :class:`~sequana.phylo.Tree` objects, all sharing the
            same set of leaf names.
        cutoff: minimum fraction of trees (0-1) a bipartition must appear in
            to be included in the consensus (0.5 = strict majority rule).

    Returns:
        Tree: consensus tree. Included internal splits carry a ``bootstrap``
        value (percentage of input trees supporting that split); leaf branch
        lengths are not meaningful in the topology-only consensus and are
        left at 0.

    Raises:
        SequanaException: if ``trees`` is empty, or the trees don't share an
            identical set of leaf names.

    Example::

        from sequana.phylo_construction import majority_consensus
        consensus = majority_consensus([tree1, tree2, tree3], cutoff=0.5)
        consensus.to_newick()
    """
    if not trees:
        raise SequanaException("majority_consensus() requires at least one tree")

    leaf_sets = [frozenset(t.leaves()) for t in trees]
    reference_leaves = leaf_sets[0]
    if any(ls != reference_leaves for ls in leaf_sets):
        raise SequanaException("All trees must share the same set of leaf names for consensus")

    n_trees = len(trees)
    split_counts: Counter = Counter()

    for tree in trees:
        for split in _get_bipartitions(tree, reference_leaves):
            split_counts[split] += 1

    accepted_splits = [(split, count / n_trees) for split, count in split_counts.items() if count / n_trees >= cutoff]
    # Larger splits (closer to the root) must be resolved before smaller
    # nested ones for the incremental tree-building below to work.
    accepted_splits.sort(key=lambda item: -len(item[0]))

    return _build_tree_from_splits(reference_leaves, accepted_splits)


def _build_tree_from_splits(all_leaves: frozenset, splits: List[Tuple[frozenset, float]]) -> Tree:
    """Incrementally build a tree from a list of (leaf-set, support) splits,
    largest splits first, by nesting each split under the smallest existing
    group that fully contains it.

    This is a simple, robust (if not maximally elegant) way to reconstruct a
    hierarchy from an accepted set of bipartitions, sufficient for
    majority-rule consensus where splits are guaranteed compatible (nested or
    disjoint) by construction.
    """
    # Each group starts as its own leaf node; groups are merged as larger
    # splits are processed.
    groups: List[Tuple[frozenset, TreeNode]] = [(frozenset([leaf]), TreeNode(name=leaf)) for leaf in sorted(all_leaves)]

    for split, support in splits:
        contained = [(leaves, node) for leaves, node in groups if leaves <= split]
        if len(contained) < 2:
            # Split already fully represented by a single existing group
            # (can happen with duplicate/nested splits); nothing to do.
            continue

        merged_leaves: frozenset = frozenset()
        for leaves, _ in contained:
            merged_leaves = merged_leaves | leaves

        if merged_leaves != split:
            # The accepted split doesn't cleanly partition the current
            # groups (can only happen with incompatible splits, which
            # majority-rule consensus above already filters out by
            # construction) -- skip defensively rather than build a wrong
            # tree.
            continue

        new_node = TreeNode(bootstrap=support * 100.0)
        for _, child in contained:
            new_node.add_child(child)

        groups = [(leaves, tnode) for leaves, tnode in groups if leaves not in [c[0] for c in contained]]
        groups.append((split, new_node))

    root = TreeNode()
    for _, node in groups:
        root.add_child(node)

    return Tree(root)
