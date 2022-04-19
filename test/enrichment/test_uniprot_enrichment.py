import json
from sequana.enrichment.uniprot_enrichment import UniprotEnrichment
import pandas as pd
import pytest

from . import test_dir


# computed
# res = pe.panther.get_enrichment(gene_lists['up'], 83333, "GO:0003674", "FISHER", "FDR")
# and saved it in json file for a mocker call
with open(f"{test_dir}/data/panther_results_up.json", "r") as fin:
    mock_values = json.loads(fin.read())


def get_gene_lists():
    up = pd.read_csv(f"{test_dir}/data/ecoli_up_gene.csv")
    down = pd.read_csv(f"{test_dir}/data/ecoli_down_gene.csv")
    up = list(up.Name)
    down = list(down.Name)
    up = up[0:200]
    down = down[0:200]
    return {'up': up, 'down': down, 'all': up+down}


def test_uniprot(mocker, tmpdir):
    gene_lists = get_gene_lists()

    pe = UniprotEnrichment(gene_lists, fc_threshold=1, taxon=83333)

    pe.compute_enrichment(ontologies=["MF"])

    df = pe.plot_go_terms("up", ontologies='MF', compute_levels=False)
    df = pe.plot_go_terms("up", ontologies='MF', compute_levels=False, log=True, show_pvalues=True)
    df = pe.plot_go_terms("up", ontologies='MF', compute_levels=False, log=True, include_negative_enrichment=True)

    pe.plot_piechart(df)

    outpng = tmpdir.join('test.png')
    pe.save_chart(df, outpng)
