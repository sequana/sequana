import os

import pytest

try:
    # Disable tqdm's background monitor thread to avoid a segfault in Python 3.13
    # caused by a race condition between the monitor thread and the GC during
    # pytest failure reporting (ast.parse triggers GC while monitor thread runs).
    import tqdm

    tqdm.tqdm.monitor_interval = 0
except ImportError:
    pass

skiptravis = pytest.mark.skipif("TRAVIS_PYTHON_VERSION" in os.environ, reason="On travis")


def pytest_runtest_setup(item):
    if "TRAVIS_PYTHON_VERSION" in os.environ:
        print("downloading toydb data from github")
        from sequana import sequana_config_path
        from sequana.taxonomy import NCBITaxonomy

        # HOME = os.getenv('HOME')
        n = NCBITaxonomy(
            "https://raw.githubusercontent.com/sequana/data/main/kraken_toydb/taxonomy/names.dmp",
            "https://raw.githubusercontent.com/sequana/data/main/kraken_toydb/taxonomy/nodes.dmp",
        )
        n.create_taxonomy_file(sequana_config_path + os.sep + "taxonomy.csv.gz")
