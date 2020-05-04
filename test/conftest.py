import os
import pytest


skiptravis = pytest.mark.skipif("TRAVIS_PYTHON_VERSION" in os.environ,
    reason="On travis")


def pytest_runtest_setup(item):
    print("hello")
    if "TRAVIS_PYTHON_VERSION" in os.environ or 1 == 1:
        from sequana.taxonomy import NCBITaxonomy
        from sequana import sequana_config_path
        n = NCBITaxonomy(
            "https://raw.githubusercontent.com/sequana/data/master/kraken_toydb/taxonomy/names.dmp",
            "https://raw.githubusercontent.com/sequana/data/master/kraken_toydb/taxonomy/nodes.dmp")
        n.create_taxonomy_file(sequana_config_path + os.sep + "taxonomy_test.dat")

