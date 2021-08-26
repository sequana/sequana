from sequana import PantherEnrichment
from easydev import TempFile
import pandas as pd
import pytest

from . import test_dir

@pytest.mark.xfail(reason="too slow or service may be down")
def test_panther():
    up = pd.read_csv(f"{test_dir}/data/ecoli_up_gene.csv")
    down = pd.read_csv(f"{test_dir}data/ecoli_down_gene.csv")  
    up = list(up.Name)
    down = list(down.Name)

    up = up[0:200]
    down = down[0:200]

    gene_lists = {'up': up, 'down': down, 'all': up+down}

    pe = PantherEnrichment(gene_lists, log2_fc_threshold=0, taxon=83333)
    pe = PantherEnrichment(gene_lists, fc_threshold=1, taxon=83333)

    pe.compute_enrichment(ontologies=["GO:0003674"],
        correction="bonferroni")

    # too slow or fails ith uniprot 404 from time to time
    #pe.get_functional_classification(pe.mygenes_up, 83333)

    ontologies = ["GO:0003674", "GO:0008150", "GO:0005575"]

    pe.compute_enrichment(  ontologies=ontologies)

    df = pe.plot_go_terms("up", ontologies=ontologies, compute_levels=False)
    df = pe.plot_go_terms("up", ontologies=ontologies, compute_levels=False, log=True)
    df = pe.plot_go_terms("up", ontologies=ontologies, compute_levels=False, log=True,
        include_negative_enrichment=True)

    pe.plot_piechart(df) 

