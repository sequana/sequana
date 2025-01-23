import pytest

from sequana.rnafold import *


def test_rnafold():
    assert pytest.approx(calculate_mfe("ACCGGTTGGTGGTGTGTGAG")) == -2
