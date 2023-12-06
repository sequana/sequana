import pandas as pd
import pytest

from sequana.viz.pca import PCA

from . import test_dir


@pytest.mark.timeout(10)
def test_pca():

    data = f"{test_dir}/data/test_pca.csv"
    df = pd.read_csv(data)
    df = df.set_index("Id")
    p = PCA(df, colors={"A1": "r", "A2": "r", "A3": "r", "B1": "b", "B2": "b", "B3": "b"})

    p.plot(n_components=2, switch_y=True, adjust=False)
    p.plot(n_components=2, switch_x=True, adjust=False)
    p.plot(n_components=3, switch_z=True, adjust=False)
    p.plot_pca_vs_max_features(n_components=4)
    p.plot_pca_vs_max_features(step=50000)
