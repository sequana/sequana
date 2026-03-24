import pandas as pd
import pytest

from sequana.viz.ideogram import Ideogram


def test_ideogram_basic_plot():
    df = pd.DataFrame(
        {
            "centromere": [5e6, 6e6],
            "length": [20e6, 22e6],
            "LHS_telomere": [2e5, 3e5],
            "RHS_telomere": [19e6, 21e6],
        }
    )
    ig = Ideogram(df)
    ig.plot()  # must not raise


def test_ideogram_with_gaps():
    df = pd.DataFrame(
        {
            "centromere": [5e6],
            "length": [20e6],
            "LHS_telomere": [2e5],
            "RHS_telomere": [19e6],
        }
    )
    ig = Ideogram(df, gaps={"1": [None, 8e6, 12e6]})
    ig.plot()
