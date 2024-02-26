import gzip

from sequana.taxonomy import NCBITaxonomy, Taxonomy

from . import test_dir


def test_taxonomy(tmp_path):
    n = NCBITaxonomy(f"{test_dir}/data/names_filtered.dmp", f"{test_dir}/data/nodes_filtered.dmp")
    filename = tmp_path / "taxo.csv.gz"
    n.create_taxonomy_file(filename)

    # but we can define another one:
    tax = Taxonomy(filename)
    tax.load_records()

    assert tax.get_lineage(2732408) == ["root", "Viruses", "Riboviria", "Orthornavirae", "Pisuviricota"]
    assert tax.get_lineage(-10) == ["unknown_taxon:-10"]
    assert tax.get_parent_taxon(2732408) == 2732396
    tax.get_lineage_and_rank(2732408)
    assert tax.get_children(2732408) == [2732506]
    assert tax.get_parent_name(2732408) == "Orthornavirae"
    assert tax.get_ranks().any()
    assert len(tax.get_records_for_given_rank("family")) == 2

    # FIXME this test works in python shell not in pytest
    assert len(tax) == 22
    # test the getter
    tax[11234]

    assert (tax.get_names_for_given_rank("family") == ["Coronaviridae", "Paramyxoviridae"]).all()

    tax.fetch_by_name("corona")

    tax.find_taxon("10684")

    # test wrong extesnion
    filename = tmp_path / "taxo.csv"
    try:
        n.create_taxonomy_file(filename)
        assert False
    except ValueError:
        assert True

    ret = tax.fetch_by_id("10090")
    ret["name"]
