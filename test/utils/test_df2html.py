import pandas as pd

from sequana.utils.df2html import df2html


def test_df2html():
    df = pd.DataFrame([[1, 2], [1, 2]])
    df2html(df)
