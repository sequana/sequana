from sequana.viz import Volcano
import pandas as pd


def test1():
    import numpy as np

    fc = np.random.randn(1000)
    pvalue = np.random.randn(1000)
    df = pd.DataFrame({"log2FoldChange": fc, "padj": pvalue ** 2})
    v = Volcano(data=df, pvalue_threshold=3)
    v.plot()
    v.plot(logy=True)
