

from sequana.utils.df2html  import df2html
import pandas as pd


def test_df2html():
    df = pd.DataFrame([[1,2], [1,2]])
    df2html(df)



