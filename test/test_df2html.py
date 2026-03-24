import pandas as pd

from sequana.utils.df2html import df2html


def test_df2html_basic():
    df = pd.DataFrame({"A": [1, 2], "B": [3, 4]})
    html = df2html(df, name="test_table")
    assert "test_table" in html
    assert "<table" in html


def test_df2html_with_links():
    df = pd.DataFrame({"name": ["a", "b"], "name_links": ["http://a.com", "http://b.com"], "other": [1, 2]})
    # df2html should identify name_links and link them to name
    html = df2html(df)
    assert "link" in html.lower()


def test_df2html_no_name():
    df = pd.DataFrame({"A": [1]})
    html = df2html(df)
    assert "<table" in html
