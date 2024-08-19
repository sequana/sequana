import pandas as pd

from sequana.viz.plotly import BinaryPercentage


def test_BinaryPercentage():
    hb = BinaryPercentage()
    hb.df = pd.DataFrame({"A": [1, 50, 100], "B": [1, 50, 100]})
    hb.plot_horizontal_bar(html_code=True)
