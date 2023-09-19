from sequana.enrichment.kegg import KEGGPathwayEnrichment
import pandas as pd
import pytest
import os
import glob
from . import test_dir



def test_ke(tmpdir):
    up = pd.read_csv(f"{test_dir}/data/ecoli_up_gene.csv")
    down = pd.read_csv(f"{test_dir}/data/ecoli_down_gene.csv")
    up = list(up.Name)
    down = list(down.Name)
    gene_lists = {'up': up, 'down': down, 'all': up +down}

    from sequana import logger
    logger.setLevel('INFO')
    ke = KEGGPathwayEnrichment(gene_lists, "eco",
            preload_directory=f"{test_dir}/data/kegg_pathways/")

    with pytest.raises(ValueError):
        ke.barplot('dummy')

    ke.barplot('up')
    ke.barplot('down')
    ke.plot_genesets_hist()
    ke.scatterplot('down')
    assert ke.find_pathways_by_gene("moaA")


    # save one pathway (just one to speed up things)
    outpng = tmpdir.join('test.png')
    df = pd.read_csv(f"{test_dir}/data/ecoli_all_gene.csv", index_col=0)
    ke.save_pathway('eco04122', df, filename=outpng)

    # save all pathways (same as input)
    path = tmpdir.mkdir("pathways_tmp")
    ke.export_pathways_to_json(str(path))

