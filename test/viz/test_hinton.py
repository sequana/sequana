def test1():
    import numpy as np

    from sequana.viz.hinton import hinton

    df = np.random.rand(20, 20) - 0.5
    hinton(df)
