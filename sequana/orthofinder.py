import glob
import itertools
from pathlib import Path

import pandas as pd
from pylab import xlabel, ylabel
from tqdm import tqdm

from sequana.lazy import numpy as np


class OrthoGroups:
    def __init__(self, directory):
        self.directory = Path(directory)
        self.ortho_groups = pd.read_csv(self.directory / "Orthogroups.tsv", sep="\t")
        self.gene_counts = pd.read_csv(self.directory / "Orthogroups.GeneCount.tsv", sep="\t")
        self.gene_counts.columns = [x + "_size" for x in self.gene_counts.columns]
        self.single_copy_orthologs = pd.read_csv(
            self.directory / "Orthogroups_SingleCopyOrthologues.txt", sep="\t", header=None
        )
        self.df_unassigned = pd.read_csv(self.directory / "Orthogroups_UnassignedGenes.tsv", sep="\t")

        # combined orthogroups and gene counts
        self.df = pd.concat([self.ortho_groups, self.gene_counts], axis=1)
        del self.gene_counts
        del self.ortho_groups

    def summary(self):
        N = len(self.df)
        print(f"The orthogroup contains {N} groups. Each sample has:")
        for x in self.df.columns:
            if x == "Orthogroup" or x.endswith("_size"):
                continue

            n = len(self.df[x].dropna())
            N = self.df[f"{x}_size"].sum()
            print(f"Found {n} groups in sample {x} in including {N} genes")

        print(f"\nUnassigned genes.")
        for sample in self.df_unassigned.columns:
            N = len(self.df_unassigned[sample].dropna())
            print(f"Found {N} unassigned genes in sample {sample}")


class GeneTrees:
    def __init__(self, directory):
        """

        Compute mean distances and variances within each group.

        Large variance with moderate mean often indicates:
            - Ancient duplication + recent paralogs
            - Asymmetric evolution across species
        """

        self.directory = Path(directory)

    def _run(self):
        from ete3 import Tree

        names = list(self.directory.glob("*tree.txt"))

        sizes = []
        mean_distances = []
        variances = []
        for name in names:
            t = Tree(str(name), format=1)
            leaves = t.get_leaves()
            n = len(leaves)
            sizes.append(n)

            if n < 2:
                mean_distances.append(0.0)
                variances.append(0.0)
                continue

            pairwise_distances = [t.get_distance(a, b) for a, b in itertools.combinations(leaves, 2)]

            mean_distances.append(np.mean(pairwise_distances))
            variances.append(np.var(pairwise_distances))

        df = pd.DataFrame(
            {"variances": variances, "sizes": sizes, "mean_distances": mean_distances, "names": [x.name for x in names]}
        )
        return df


class OrthoFinder:
    """

    o = OrthoFinder(".")
    g = GFF3("Ld1S.gff")
    annotations, chromosomes = o.add_annotation_Ld1S(g)
    o.ortho_groups['annotation'] = annotations
    o.ortho_groups['chromosome'] = chromosomes


    """

    def __init__(self, directory):
        self.directory = Path(directory)
        # this contains all orthogroup names for each species.
        self.ortho_groups = pd.read_csv(self.directory / "Orthogroups" / "Orthogroups.tsv", sep="\t")
        self.gene_counts = pd.read_csv(self.directory / "Orthogroups" / "Orthogroups.GeneCount.tsv", sep="\t")

    def summary(self):
        N = len(self.ortho_groups)
        print(f"The orthogroup contains {N} groups")
        for x in self.ortho_groups.columns:
            if x == "Orthogroup":
                continue
            print(x, len(self.ortho_groups[x].dropna()))

    def hist_gene_counts(self, bins=50):
        self.gene_counts.Total.hist(bins=bins, log=True)
        xlabel("Total gene count per group")
        ylabel("#")

    def add_annotation_Ld1S(self, gff, column="Ld1S", genetic_type="gene"):

        # Build lookup table once
        df = gff.df.query("genetic_type == @genetic_type").copy()
        df = df[["gene_id", "combinedAnnotation", "seqid", "start"]]
        lookup = df.set_index("gene_id").to_dict("index")

        annotations = []
        chromosomes = []
        positions = []

        for IDs in tqdm(self.ortho_groups["Ld1S"]):

            try:
                genes = [x.strip(",") for x in IDs.split()]

                hits = [lookup[g] for g in genes if g in lookup]

                if hits:
                    annotations.append(hits[0]["combinedAnnotation"])
                    chromosomes.append(sorted(list(set([h["seqid"] for h in hits]))))
                    positions.append(" ".join([str(h["start"]) for h in hits]))
                else:
                    annotations.append("")
                    chromosomes.append("")
                    positions.append("")

            except Exception:
                annotations.append("")
                chromosomes.append("")
                positions.append("")
        return annotations, chromosomes, positions
