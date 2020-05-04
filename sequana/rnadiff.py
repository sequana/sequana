from pathlib import Path
import re
import os

from sequana.lazy import pandas as pd
from sequana.lazy import pylab

from matplotlib_venn import venn2_unweighted, venn3_unweighted
try:
    import gseapy
except:
    pass


class RNADiffResults(object):
    """ A representation of a biomics RNADiffResults analysis.
    """

    TABLE_DIR_NAME = "tables"
    # File name for the complete table is of form:
    # B1234-v1 (Number of project, number of version)
    TOTAL_TABLE_PATTERN = ".*B\d+-v\d+\.complete\.xls"
    NORMCOUNTS_TABLE = ".*B\d+-v\d+\.normCounts\.xls"
    SEP = "\t"

    def __init__(self, rnadiff_folder, alpha=0.05, out_dir="gsea"):
        self.analysis_folder = Path(rnadiff_folder)
        self.name = self.analysis_folder.stem
        self.table_folder = self.analysis_folder / self.TABLE_DIR_NAME
        self.df = self.get_table(self.TOTAL_TABLE_PATTERN)
        self.normcounts = self.get_table(self.NORMCOUNTS_TABLE)
        self.comparisons = self.get_comparisons()
        self.dr_gene_lists = self.get_gene_lists(alpha=alpha)
        self.out_dir = out_dir
        self.alpha = alpha

    def get_table(self, pattern):
        """Extract the complete (with all comparisons) table from a biomics
        RNADiff analysis or the normCounts table depending on the PATTERN specified.
        """

        table_files = [f for f in self.table_folder.glob("*.xls")]

        table = [f for f in table_files if re.match(pattern, str(f))]

        if len(table) == 1 and table[0].is_file():

            return pd.read_csv(table[0], self.SEP, index_col=0)

        else:
            raise IOError(
                f"Cannot find a single proper table with pattern: {pattern} from RNADiff: {table}"
            )

    def get_comparisons(self):

        comparisons = set(
            [x.split(".")[0] for x in self.df.filter(regex=(".*vs.*\..*")).columns]
        )
        return comparisons

    def get_gene_lists(self, alpha=None):

        gene_sets = {}

        for compa in self.comparisons:
            gene_sets[compa] = {}
            regex = f"{compa}\..*"

            log2FC_colname = f"{compa}.log2FoldChange"
            padj_colname = f"{compa}.padj"
            sub_df = self.df.filter(regex=regex)

            query_up = f"`{log2FC_colname}` > 0 and `{padj_colname}` < {alpha}"
            query_down = f"`{log2FC_colname}` < 0 and `{padj_colname}` < {alpha}"
            query_all = f"`{padj_colname}` < {alpha}"

            gene_sets[compa]["up"] = set(sub_df.query(query_up).index)
            gene_sets[compa]["down"] = set(sub_df.query(query_down).index)
            gene_sets[compa]["all"] = set(sub_df.query(query_all).index)

        return gene_sets

    def run_enrichr(self, gene_sets, top_term=20):
        """
        Run enrichr using gseapy

        gene_sets: A list of databases coming from the database available with
        gseapy (see gseapy.get_library_name(database="Mouse"))
        top_term: The max number of term to plot.
        """

        out_dir_enr = os.path.join(self.out_dir, "enrichr")
        os.makedirs(out_dir_enr, exist_ok=True)
        failed_comparisons = []
        no_dr_genes = []

        for compa in self.dr_gene_lists:
            for direction in self.dr_gene_lists[compa]:

                ensembl_ids = self.dr_gene_lists[compa][direction]
                gene_symbols = [x.upper() for x in self.df.loc[ensembl_ids].gene_name]

                if len(gene_symbols) == 0:
                    no_dr_genes.append(f"{compa}_{direction}")
                    continue

                for gene_set in gene_sets:

                    # print(f"{gene_set} {compa} {direction} {len(gene_symbols)} ")

                    enr = gseapy.enrichr(
                        gene_list=gene_symbols, gene_sets=gene_set, no_plot=True
                    )

                    title = f"{gene_set}_{compa}_{direction}"

                    out_dir_db = os.path.join(out_dir_enr, gene_set)
                    os.makedirs(out_dir_db, exist_ok=True)

                    # gseapy.plot.barplot(enr, cutoff=0.05)
                    try:
                        gseapy.plot.dotplot(
                            enr.res2d,
                            cutoff=self.alpha,
                            top_term=top_term,
                            ofname=os.path.join(out_dir_db, title + ".pdf"),
                            title=title,
                        )

                    except:
                        failed_comparisons.append(f"{gene_set}_{compa}_{direction}")

        print(f"Failed for: {','.join(failed_comparisons)}")
        print(f"No DR genes found for: {','.join(no_dr_genes)}")

    def summary(self):
        summary = pd.DataFrame(
            {
                k: {x: len(self.dr_gene_lists[k][x]) for x in self.dr_gene_lists[k]}
                for k in self.dr_gene_lists
            }
        )

        return summary

    def venn(self, compa_list, direction="all", prefix=""):
        """
        Plot a venn diagram comparing the list compa_list of dr gene lists.
        compa_list is a list of comparison names from Deseq2 results
        direction specifies either if up/down/all dr genes are considered
        prefix is a string to be added as prefix to the outfile name.

        compa_list can be a list of lists of comparisons to make.
        ie [["WT", "KO1"],["WT", "KO2"]
        """

        # If compa_list is a list of lists of comparison
        if all(isinstance(l, list) for l in compa_list):

            fig, ax = pylab.subplots(6, 1, figsize=(6, 20))
            ax = ax.flat

            for i, c in enumerate(compa_list):

                plot_venn(
                    [self.dr_gene_lists[x][direction] for x in c],
                    [compa_name for compa_name in c],
                    ax=ax[i],
                )
        # If compa is only a list of comparisons
        else:
            plot_venn(
                [self.dr_gene_lists[x][direction] for x in compa_list],
                [compa_name for compa_name in compa_list],
            )
        out_dir = os.path.join(self.out_dir, "vennDiagrams")
        os.makedirs(out_dir, exist_ok=True)
        outfile = os.path.join(out_dir, f"{prefix}vennDiagrams_{direction}.pdf")

        pylab.savefig(outfile, bbox_inches="tight")

    def compare(self, rnadiff_res_obj, make_plot=False, plot_to_file=""):
        """Compare two RNADiffResults objects.  For now this will plot venn
        diagrams for the comparisons considering UP and DOWN regulated
        genes separately

        If make_plot, will plot venn diagrams

        If simple, returns the overlap (regardless of up or down regulated)
        """

        res = {}

        for compa in self.dr_gene_lists:

            # Up
            up1 = self.dr_gene_lists[compa]["up"]
            up2 = rnadiff_res_obj.dr_gene_lists[compa]["up"]

            # Down
            down1 = self.dr_gene_lists[compa]["down"]
            down2 = rnadiff_res_obj.dr_gene_lists[compa]["down"]

            res[compa] = {
                self.name + "_specific_up": self.df.loc[up1 - up2],
                rnadiff_res_obj.name
                + "_specific_up": rnadiff_res_obj.df.loc[up2 - up1],
                self.name + "_specific_down": self.df.loc[down1 - down2],
                rnadiff_res_obj.name
                + "_specific_down": rnadiff_res_obj.df.loc[down2 - down1],
                self.name
                + "_"
                + rnadiff_res_obj.name
                + "_common_up": self.df.loc[up1.intersection(up2)],
                self.name
                + "_"
                + rnadiff_res_obj.name
                + "_common_down": self.df.loc[down1.intersection(down2)],
            }

        res_summary = {}
        for compa in res:
            res_summary[compa] = {}
            for direction in res[compa]:
                res_summary[compa][direction] = res[compa][direction].shape[0]

        print(pd.DataFrame(res_summary).to_string())


def plot_venn(compa_list, labels, title=None, ax=None):
    """ Plot venn diagramm according to number of groups.
    """

    if len(compa_list) == 2:
        venn_function = venn2_unweighted
    elif len(compa_list) == 3:
        venn_function = venn3_unweighted
    else:
        raise IOError("Venn diagramm supports only 2 or 3 groups.")

    venn_function(compa_list, set_labels=labels, ax=ax)
    # ax.set_title = title
