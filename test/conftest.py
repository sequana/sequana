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

try:
    # Fix kaleido 0.2.1 segfault on Python 3.13.
    # The original _collect_standard_error does not hold a strong reference to
    # self._proc inside the loop.  When the kaleido subprocess exits, Python
    # 3.13's more-aggressive GC can collect the Popen object while readline()
    # is still executing on the now-freed file descriptor, causing a SIGSEGV.
    # Holding a local reference to `proc` keeps the object alive for the
    # duration of each readline() call and lets the loop exit cleanly on EOF.
    from kaleido.scopes.base import BaseScope

    def _safe_collect_standard_error(self):
        while True:
            proc = self._proc  # keep a strong reference to prevent GC
            if proc is None:
                return
            try:
                val = proc.stderr.readline()
            except (OSError, ValueError):
                return
            if not val:  # EOF – subprocess has terminated
                return
            self._std_error.write(val)

    BaseScope._collect_standard_error = _safe_collect_standard_error
except ImportError:
    pass

@pytest.fixture(scope="session", autouse=True)
def _shutdown_kaleido_after_session():
    """Ensure the kaleido subprocess is shut down after the test session.

    kaleido 0.2.1 keeps a Chromium subprocess alive in a module-level global
    (plotly.io.kaleido.scope).  Without an explicit shutdown the background
    stderr-reader thread keeps running until Python exits, which can cause
    a SIGSEGV on Python 3.13 due to the more-aggressive incremental GC.
    """
    yield
    try:
        from plotly.io import kaleido

        kaleido.scope._shutdown_kaleido()
    except Exception:
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
