import glob
import os

import pandas as pd
import pytest

from sequana.enrichment.kegg import KEGGPathwayEnrichment

from . import test_dir


def test_ke(tmpdir):
    up = pd.read_csv(f"{test_dir}/data/ecoli_up_gene.csv")
    down = pd.read_csv(f"{test_dir}/data/ecoli_down_gene.csv")
    up = list(up.Name)
    down = list(down.Name)
    gene_lists = {"up": up, "down": down, "all": up + down}

    # used genes should be all genes used in the analyses. We do not have it in this example,
    # so we just add up and down
    ke = KEGGPathwayEnrichment(
        gene_lists, "eco", preload_directory=f"{test_dir}/data/kegg_pathways/", used_genes=up + down
    )

    # set background ourself
    ke = KEGGPathwayEnrichment(gene_lists, "eco", preload_directory=f"{test_dir}/data/kegg_pathways/", background=10000)

    # let use use the kegg genes (no background, no used_genes)
    ke = KEGGPathwayEnrichment(
        gene_lists,
        "eco",
        preload_directory=f"{test_dir}/data/kegg_pathways/",
    )

    with pytest.raises(ValueError):
        ke.barplot("dummy")

    ke.barplot("up")
    ke.barplot("down")
    ke.plot_genesets_hist()
    ke.scatterplot("down")
    assert ke.find_pathways_by_gene("moaA")
    ke.barplot_up_and_down()

    # save one pathway (just one to speed up things)
    outpng = tmpdir.join("test.png")
    df = pd.read_csv(f"{test_dir}/data/ecoli_all_gene.csv", index_col=0)
    ke.save_pathway("eco04122", df, filename=outpng)

    # save all pathways (same as input)
    path = tmpdir.mkdir("pathways_tmp")
    ke.save_pathways(str(path))

    ke.save_project("TEST", outdir=path)
