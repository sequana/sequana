import pytest

from sequana.gfa import GFA

from . import test_dir


def test_gfa(tmp_path):
    gfa_file = tmp_path / "test.gfa"
    gfa_file.write_text(
        "S\tedge1\tACGT\n"
        "L\tedge1\t+\tedge2\t-\t0M\n"
        "H\tVN:Z:1.0\n"
    )
    # GFA just prints output; verify it runs without error
    g = GFA(str(gfa_file))
