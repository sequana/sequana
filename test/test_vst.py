import pandas as pd

from sequana.vst import VST

from . import test_dir


def test_vst():

    VST()
    data = f"{test_dir}/data/misc/test_pca.csv"
    df = pd.read_csv(data)["A1"]
    # to avoid warnings:
    df = df.replace(0, 1)

    new_df = VST.anscombe(df)
    new_df = VST.inverse_anscombe(df)
    new_df = VST.generalized_anscombe(df, 1, 1)

    df = list(df)

    new_df = VST.anscombe(df)
    new_df = VST.inverse_anscombe(df)
    new_df = VST.generalized_anscombe(df, 1, 1)
