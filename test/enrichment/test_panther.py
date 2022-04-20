import json
from sequana import PantherEnrichment
import pandas as pd
import pytest

from . import test_dir


# computed
# res = pe.panther.get_enrichment(gene_lists['up'], 83333, "GO:0003674", "FISHER", "FDR")
# and saved it in json file for a mocker call
def get_mock_data(*args, **kwargs):
    with open(f"{test_dir}/data/panther_results_up.json", "r") as fin:
        mock_values = json.loads(fin.read())
    return mock_values


def get_gene_lists():
    up = pd.read_csv(f"{test_dir}/data/ecoli_up_gene.csv")
    down = pd.read_csv(f"{test_dir}/data/ecoli_down_gene.csv")
    up = list(up.Name)
    down = list(down.Name)
    up = up[0:200]
    down = down[0:200]
    return {'up': up, 'down': down, 'all': up+down}



def test_panther1(mocker, tmpdir):
    mocker.patch('bioservices.panther.Panther.get_enrichment', side_effect=get_mock_data)
    gene_lists = get_gene_lists()

    pe = PantherEnrichment(gene_lists, log2_fc_threshold=0, taxon=83333)
    pe = PantherEnrichment(gene_lists, fc_threshold=1, taxon=83333)

    pe.compute_enrichment(ontologies=["MF"], correction='bonferroni')

    df = pe.plot_go_terms("up", ontologies='MF', compute_levels=False)
    df = pe.plot_go_terms("up", ontologies='MF', compute_levels=False, log=True, show_pvalues=True)
    df = pe.plot_go_terms("up", ontologies='MF', compute_levels=False, log=True, include_negative_enrichment=True)

    pe.plot_piechart(df)

    outpng = tmpdir.join('test.png')
    pe.save_chart(df, outpng)


def test_panther_all(mocker):
    mocker.patch('bioservices.panther.Panther.get_enrichment', side_effect=get_mock_data)

    gene_lists = get_gene_lists()

    pe = PantherEnrichment(gene_lists, fc_threshold=1, taxon=83333)

    ontologies = ["MF", "CC", "BP"]

    pe.compute_enrichment(ontologies=ontologies)

    # plots
    df = pe.plot_go_terms("up", ontologies=ontologies, compute_levels=False, log=True)

    res = pe._get_graph(['GO:0015343', 'GO:0003674', 'GO:0015344', 'GO:0051082', 'GO:0038023', 
        'GO:0016887', 'GO:0015620', 'GO:0051087', 'GO:0140657', 'GO:0005488'], 'MF')

