"""
Tree construction: UPGMA, neighbor-joining, and consensus trees
=================================================================

Build phylogenetic trees directly from a pairwise distance matrix (no
external tool required), and summarize multiple bootstrap-replicate trees
into a majority-rule consensus tree.
"""
from pylab import *

from sequana.phylo import Tree
from sequana.phylo_construction import majority_consensus, neighbor_joining, upgma

##############################################################################
# UPGMA vs neighbor-joining on the same distance matrix. UPGMA assumes a
# constant molecular clock (produces an ultrametric tree -- every leaf
# equidistant from the root); neighbor-joining does not, and is the more
# broadly applicable choice for real sequence divergence data.

names = ["human", "chimp", "gorilla", "orangutan"]
matrix = [
    [0, 2, 4, 8],
    [2, 0, 4, 8],
    [4, 4, 0, 8],
    [8, 8, 8, 0],
]

upgma_tree = upgma(matrix, names)
nj_tree = neighbor_joining(matrix, names)

print("UPGMA:", upgma_tree.to_newick())
print("NJ   :", nj_tree.to_newick())


##############################################################################
# A minimal rectangular tree plot that preserves the real topology and branch
# lengths (unlike a distance-recomputed dendrogram) -- draws directly from
# the parsed/constructed ``Tree`` object.


def plot_tree(tree, ax, title):
    leaves = tree.leaves()
    y_positions = {name: i for i, name in enumerate(leaves)}

    def node_y(node):
        if node.is_leaf():
            return y_positions[node.name]
        ys = [node_y(c) for c in node.children]
        return sum(ys) / len(ys)

    def node_x(node, depth=0.0):
        return depth + node.branch_length

    def draw(node, x_parent=0.0):
        x = x_parent + node.branch_length
        y = node_y(node)
        ax.plot([x_parent, x], [y, y], color="black")
        if not node.is_leaf():
            child_ys = [node_y(c) for c in node.children]
            ax.plot([x, x], [min(child_ys), max(child_ys)], color="black")
            for child in node.children:
                draw(child, x)
        else:
            ax.text(x + 0.15, y, node.name, va="center", fontsize=10)

    draw(tree.root)
    ax.set_yticks([])
    ax.set_xlabel("distance")
    ax.set_title(title)
    ax.set_xlim(-0.5, 10)


fig, axes = subplots(1, 2, figsize=(10, 4))
plot_tree(upgma_tree, axes[0], "UPGMA")
plot_tree(nj_tree, axes[1], "Neighbor-joining")
tight_layout()

##############################################################################
# Majority-rule consensus of bootstrap replicate trees: splits present in
# more than ``cutoff`` fraction of the input trees are retained, annotated
# with their support percentage.

t1 = Tree.from_newick("((human,chimp),(gorilla,orangutan));")
t2 = Tree.from_newick("((human,chimp),(gorilla,orangutan));")
t3 = Tree.from_newick("((human,gorilla),(chimp,orangutan));")

consensus = majority_consensus([t1, t2, t3], cutoff=0.5)
print("Consensus:", consensus.to_newick())
print(consensus.to_ascii())
