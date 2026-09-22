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


@pytest.fixture(autouse=True)
def reset_config_output_dir():
    """Reset config.output_dir and clean up css/js/images dirs created by tests."""
    import shutil

    from sequana.utils import config

    original_output_dir = config.output_dir
    yield
    config.output_dir = original_output_dir

    # Clean up css/js/images directories if they were created in the current working directory
    cwd = os.getcwd()
    for dirname in ("css", "js", "images"):
        path = os.path.join(cwd, dirname)
        if os.path.isdir(path):
            try:
                shutil.rmtree(path)
            except OSError:
                # Handle case where directory is in use or contains locked files
                pass


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
