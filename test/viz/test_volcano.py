from sequana.viz import Volcano

def test1():
    import numpy as np
    fc = np.random.randn(1000)
    pvalue = np.random.randn(1000)
    v = Volcano(fc, -np.log10(pvalue**2), pvalue_threshold=3)
    v.plot()
    v.plot(logy=True)


