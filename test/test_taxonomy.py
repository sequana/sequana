from sequana.taxonomy import Taxonomy
from sequana.taxonomy import NCBITaxonomy

from . import test_dir

def test_taxonomy(tmp_path):
    n = NCBITaxonomy(f"{test_dir}/data/names_filtered.dmp", f"{test_dir}/data/nodes_filtered.dmp")
    filename = tmp_path / "taxo.dat"
    n.create_taxonomy_file(filename)

    # default is to load local file in config/sequana
    #tax = Taxonomy()

    # but we can define another one:
    tax = Taxonomy(filename)
    tax.records = {}
    tax.load_records()
    tax.load_records()


    assert  tax.get_lineage(2732408)==['root', 'Viruses', 'Riboviria', 'Orthornavirae', 'Pisuviricota']
    assert  tax.get_lineage(-10) == ['unknown_taxon:-10']
    assert tax.get_parent_taxon(2732408) == 2732396
    tax.get_lineage_and_rank(2732408)
    assert  tax.get_children(2732408) == []
    assert  tax.get_parent_name(2732408) == "Orthornavirae"
    assert tax.get_ranks()
    assert len(tax.get_record_for_given_rank("family")) == 2
    

    #FIXME this test works in python shell not in pytest
    assert len(tax) == 22
    # test the getter
    tax[11234]


    assert tax.get_names_for_given_rank("family") ==  ['Coronaviridae', 'Paramyxoviridae']

    tax.append_existing_database(filename)

    tax.fetch_by_name('corona')

    tax.find_taxon('10684')
