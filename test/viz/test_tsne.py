import pandas as pd
import pytest

from sequana.viz.tsne import TSNE

from . import test_dir


@pytest.mark.timeout(10)
def test_tsne():

    data = f"{test_dir}/data/test_pca.csv"
    df = pd.read_csv(data)
    df = df.set_index("Id")
    p = TSNE(df, colors={"A1": "r", "A2": "r", "A3": "r", "B1": "b", "B2": "b", "B3": "b"})

    p.plot(perplexity=2)
