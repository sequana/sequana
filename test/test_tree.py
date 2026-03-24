from pathlib import Path

import pytest

from sequana.utils.tree import DisplayablePath, HTMLDirectory


def test_displayable_path(tmp_path):
    # Create a small dummy tree
    (tmp_path / "dir1").mkdir()
    (tmp_path / "dir1" / "file1.txt").touch()
    (tmp_path / "file2.txt").touch()

    paths = list(DisplayablePath.make_tree(tmp_path))
    assert len(paths) >= 3
    for p in paths:
        assert isinstance(p.displayable(), str)


def test_html_directory_no_filter(tmp_path):
    (tmp_path / "test.html").touch()
    hd = HTMLDirectory(str(tmp_path))
    html = hd.get_html()
    assert "<html>" in html
    assert "test.html" in html


def test_html_directory_with_pattern(tmp_path):
    (tmp_path / "match.txt").touch()
    (tmp_path / "nomatch.txt").touch()
    hd = HTMLDirectory(str(tmp_path), pattern="match.txt")
    html = hd.get_html()
    assert "match.txt" in html
    assert "nomatch.txt" not in html


def test_html_directory_with_skip(tmp_path):
    (tmp_path / "skip_me.txt").touch()
    (tmp_path / "keep_me.txt").touch()
    hd = HTMLDirectory(str(tmp_path), skip_pattern=["skip_me"])
    html = hd.get_html()
    assert "keep_me.txt" in html
    assert "skip_me.txt" not in html
