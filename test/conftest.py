import os
import pytest


skiptravis = pytest.mark.skipif("TRAVIS_PYTHON_VERSION" in os.environ,
    reason="On travis")


def pytest_runtest_setup(item):
    if "TRAVIS_PYTHON_VERSION" in os.environ:
        print("downloading toydb data from github")
        from sequana.taxonomy import NCBITaxonomy
        from sequana import sequana_config_path
        #HOME = os.getenv('HOME')
        n = NCBITaxonomy(
            "https://raw.githubusercontent.com/sequana/data/main/kraken_toydb/taxonomy/names.dmp",
            "https://raw.githubusercontent.com/sequana/data/main/kraken_toydb/taxonomy/nodes.dmp")
        n.create_taxonomy_file(sequana_config_path + os.sep + "taxonomy.dat")

