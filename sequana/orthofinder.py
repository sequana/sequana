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

        print("\nUnassigned genes.")
        for sample in self.df_unassigned.columns:
            if sample == "Orthogroup":
                continue
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
    annotations, chromosomes, positions = o.add_annotation_Ld1S(g)
    o.ortho_groups['annotation'] = annotations
    o.ortho_groups['chromosome'] = chromosomes
    o.ortho_groups['position'] = positions


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

        for IDs in tqdm(self.ortho_groups[column]):
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

            except AttributeError:
                # IDs is NaN (empty cell in the orthogroups table for this
                # species) -- pandas represents it as float('nan'), which has
                # no .split() method.
                annotations.append("")
                chromosomes.append("")
                positions.append("")
        return annotations, chromosomes, positions


class OrthoFinderResults:
    """Summaries of an OrthoFinder results directory (data only, no plotting).

    The directory must contain ``Orthogroups/`` and may contain
    ``Comparative_Genomics_Statistics/`` and ``Species_Tree/`` (as produced by
    OrthoFinder). Used by the sequana_orthofinder pipeline for its report.

    ::

        res = OrthoFinderResults("ResultsLight")
        res.pairwise_shared()
        res.intersections()
        res.category_counts()
    """

    RELATIONS = ["one-to-one", "one-to-many", "many-to-one", "many-to-many"]

    def __init__(self, directory):
        self.directory = Path(directory)
        og = self.directory / "Orthogroups"
        self.orthogroups = pd.read_csv(og / "Orthogroups.tsv", sep="\t", index_col=0)
        self.gene_counts = pd.read_csv(og / "Orthogroups.GeneCount.tsv", sep="\t", index_col=0)
        self.species = list(self.orthogroups.columns)
        self.presence = self.gene_counts[self.species] > 0

    # ------------------------------------------------------------------ gene content
    def pairwise_shared(self):
        """Orthogroups shared by each pair of genomes.

        Percentage is relative to orthogroups present in either genome. This
        is gene-content sharing, not gene-order (synteny) conservation.
        """
        rows = []
        for a, b in itertools.combinations(self.species, 2):
            shared = int((self.presence[a] & self.presence[b]).sum())
            union = int((self.presence[a] | self.presence[b]).sum())
            rows.append((a, b, shared, union, round(100 * shared / union, 2) if union else 0.0))
        return pd.DataFrame(rows, columns=["Genome1", "Genome2", "Shared_Orthogroups", "Union_Orthogroups", "Shared_%"])

    def intersections(self):
        """Number of orthogroups for each combination of genomes (UpSet data).

        :return: DataFrame with one boolean column per genome plus ``count``,
            sorted by decreasing count.
        """
        df = self.presence.groupby(self.species).size().rename("count").reset_index()
        return df.sort_values("count", ascending=False).reset_index(drop=True)

    def categories(self):
        """Classify each orthogroup as core, shell or specific.

        core: present in all genomes; specific: present in a single genome;
        shell: anything in between (with a single genome everything is core).
        """
        n = self.presence.sum(axis=1)
        cat = pd.Series("shell", index=self.presence.index, name="category")
        cat[n == len(self.species)] = "core"
        if len(self.species) > 1:
            cat[n == 1] = "specific"
        return cat

    def category_counts(self):
        """Per-genome count of orthogroups in each category (DataFrame genome x category)."""
        cat = self.categories()
        out = (
            pd.DataFrame({sp: cat[self.presence[sp]].value_counts() for sp in self.species})
            .T.reindex(columns=["core", "shell", "specific"])
            .fillna(0)
            .astype(int)
        )
        out.index.name = "Genome"
        return out

    def copy_number_distribution(self, max_copies=10):
        """Orthogroup count per gene copy number, per genome (columns), capped at ``max_copies``."""
        capped = self.gene_counts[self.species].clip(upper=max_copies)
        out = pd.DataFrame({sp: capped[sp][capped[sp] > 0].value_counts() for sp in self.species})
        return out.sort_index().fillna(0).astype(int)

    def expanded_families(self, top=None, with_genes=False):
        """Orthogroups ranked by copy-number spread across genomes.

        Only orthogroups present in all genomes are considered, so that the
        spread reflects expansion/contraction and not presence/absence.

        :param with_genes: add the columns ``Largest_in`` (genome with the most
            copies) and ``Genes`` (its gene IDs, comma-separated)
        """
        counts = self.gene_counts[self.species]
        df = counts[self.presence.all(axis=1)].copy()
        df["Min"] = df.min(axis=1)
        df["Max"] = df.max(axis=1)
        df["Range"] = df["Max"] - df["Min"]
        df["Variance"] = counts.loc[df.index].var(axis=1).round(2)
        df = df.sort_values(["Range", "Variance"], ascending=False)
        df.index.name = "Orthogroup"
        df = df.head(top) if top else df
        if with_genes:
            df["Largest_in"] = counts.loc[df.index].idxmax(axis=1)
            df["Genes"] = [self.orthogroups.loc[og, sp] for og, sp in zip(df.index, df["Largest_in"])]
        return df

    def specific_orthogroups(self):
        """Orthogroups present in a single genome: Orthogroup, Genome, N_genes, Genes."""
        cat = self.categories()
        rows = []
        if len(self.species) > 1:
            for og in cat.index[cat == "specific"]:
                sp = self.presence.columns[self.presence.loc[og].values.argmax()]
                rows.append((og, sp, int(self.gene_counts.loc[og, sp]), self.orthogroups.loc[og, sp]))
        return pd.DataFrame(rows, columns=["Orthogroup", "Genome", "N_genes", "Genes"])

    # ------------------------------------------------------------------ OrthoFinder statistics
    @property
    def statistics_dir(self):
        return self.directory / "Comparative_Genomics_Statistics"

    def _read(self, name, **kwargs):
        path = self.statistics_dir / name
        return pd.read_csv(path, sep="\t", **kwargs) if path.exists() else None

    def statistics_overall(self):
        """Statistics_Overall.tsv as a two-column (Statistic, Value) table of the numeric part."""
        path = self.statistics_dir / "Statistics_Overall.tsv"
        if not path.exists():
            return None
        rows = []
        for line in path.read_text().splitlines():
            if not line.strip():
                break  # first block only; the following ones are histograms
            key, _, value = line.partition("\t")
            rows.append((key, value))
        return pd.DataFrame(rows, columns=["Statistic", "Value"])

    def statistics_per_species(self):
        """Statistics_PerSpecies.tsv as a (Statistic x genome) table of the numeric part."""
        path = self.statistics_dir / "Statistics_PerSpecies.tsv"
        if not path.exists():
            return None
        lines = []
        for line in path.read_text().splitlines():
            if not line.strip():
                break
            lines.append(line.split("\t"))
        df = pd.DataFrame(lines[1:], columns=["Statistic"] + lines[0][1:])
        return df

    def duplications_per_node(self):
        """Gene duplications mapped to each species-tree node (None if unavailable)."""
        return self._read("Duplications_per_Species_Tree_Node.tsv")

    def ortholog_relations(self):
        """Orthology relations for each ordered genome pair (A -> B).

        Columns: Genome1, Genome2, one-to-one, one-to-many, many-to-one,
        many-to-many, Total (number of orthologues, as counted by OrthoFinder).
        """
        mats = {rel: self._read(f"OrthologuesStats_{rel}.tsv", index_col=0) for rel in self.RELATIONS}
        mats["Total"] = self._read("OrthologuesStats_Totals.tsv", index_col=0)
        if any(m is None for m in mats.values()):
            return None
        rows = []
        for a in self.species:
            for b in self.species:
                if a != b:
                    rows.append((a, b) + tuple(int(mats[k].loc[a, b]) for k in self.RELATIONS + ["Total"]))
        return pd.DataFrame(rows, columns=["Genome1", "Genome2"] + self.RELATIONS + ["Total"])

    # ------------------------------------------------------------------ species tree
    @property
    def species_tree_file(self):
        """OrthoFinder rooted species tree (None for fewer than 3 genomes)."""
        path = self.directory / "Species_Tree" / "SpeciesTree_rooted_node_labels.txt"
        return path if path.exists() else None


def gene_tree_discordance(gene_trees_dir, species, species_tree):
    """Normalised Robinson-Foulds distance between gene trees and the species tree.

    Only single-copy orthogroups (exactly one gene per genome) can be compared
    directly; the others are skipped. At least 4 genomes are needed for the
    distance to be informative.

    :param gene_trees_dir: OrthoFinder ``Gene_Trees`` directory
    :param species: genome names; gene tree leaves are ``<genome>_<gene id>``
    :param species_tree: Newick file of the species tree
    :return: Series indexed by orthogroup name (``OG0000000``) with values in [0, 1]
    """
    from ete3 import Tree

    if len(species) < 4:
        return pd.Series(dtype=float, name="RF")

    prefixes = sorted(species, key=len, reverse=True)  # genome names may contain "_"

    def to_species(leaf):
        for sp in prefixes:
            if leaf.startswith(sp + "_"):
                return sp
        return None

    reference = Tree(str(species_tree), format=1)
    for node in reference.traverse():
        if not node.is_leaf():
            node.name = ""
    values = {}
    for path in Path(gene_trees_dir).glob("*_tree.txt"):
        tree = Tree(str(path), format=1)
        leaves = tree.get_leaves()
        names = [to_species(leaf.name) for leaf in leaves]
        if len(leaves) != len(species) or sorted(map(str, names)) != sorted(species):
            continue
        for leaf, name in zip(leaves, names):
            leaf.name = name
        rf, max_rf = tree.robinson_foulds(reference.copy(), unrooted_trees=True)[:2]
        if max_rf:
            values[path.name.replace("_tree.txt", "")] = rf / max_rf
    return pd.Series(values, name="RF", dtype=float)
