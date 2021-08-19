from sequana.enrichment.kegg import KEGGPathwayEnrichment
import pandas as pd
import pytest
import os
import glob
from . import test_dir



@pytest.mark.xfail(reason="too slow or service may be down")
def test_ke():
    up = pd.read_csv(f"{test_dir}/data/ecoli_up_gene.csv")
    down = pd.read_csv(f"{test_dir}/data/ecoli_down_gene.csv")
    up = list(up.Name)
    down = list(down.Name)
    gene_lists = {'up': up, 'down': down, 'all': up +down}

    ke = KEGGPathwayEnrichment(gene_lists, "lbi", log2_fc=0)
    ke.barplot('down')
    ke.barplot('up')
    ke.barplot('down')
    ke.plot_genesets_hist()
    ke.scatterplot('down')
    assert ke.find_pathways_by_gene("moaA")
    assert ke.find_pathways_by_gene("moaA", match="exact")

    # cleanup
    filenames = glob.glob(f"{test_dir}/test/CUSTOM*")
    for filename in filenames:
        os.remove(filename)
    os.removedirs(f"{test_dir}/test")
 

