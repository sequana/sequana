import pandas as pd

from sequana.viz import Volcano


def test1():
    import numpy as np

    fc = np.random.randn(50)
    pvalue = np.random.randn(50)
    df = pd.DataFrame(
        {"log2FoldChange": fc, "padj": (pvalue**2) / 100.0, "annot": ["gene" + str(i) for i in range(50)]}
    )
    v = Volcano(data=df, pvalue_threshold=1)
    v.plot()
    v.plot(logy=True)
    v.annotate()
