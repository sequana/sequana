from sequana.taxonomy import Taxonomy


def test_taxonomy():
    from sequana.taxonomy import NCBITaxonomy                                                        
    n = NCBITaxonomy("https://raw.githubusercontent.com/sequana/data/master/kraken_toydb/taxonomy/names.dmp", "https://raw.githubusercontent.com/sequana/data/master/kraken_toydb/taxonomy/nodes.dmp")
    n.create_taxonomy_file("taxo.dat")           
