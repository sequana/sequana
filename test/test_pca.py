from sequana.viz.pca import PCA


def test_pca():

    from sequana import sequana_data 
    import pandas as pd
    filename = sequana_data("test_pca.csv")
    df = pd.read_csv(filename)
    df = df.set_index("Id")

    p = PCA(df)
    p = PCA(df, colors={
        "A1": "r", "A2": "r", "A3": "r",
        "B1": "b", "B2": "b", "B3": "b"})

    # default
    p.plot()

    # test some optins
    p.plot(transform="log", switch_y=True, switch_x=True, max_features=500)
    p.plot(transform="anscombe", switch_y=True, switch_x=True, max_features=500)

    # set max_features to large values to test option
    p.plot_pca_vs_max_features(step=10000)
    p.plot_pca_vs_max_features(n_components=3)

    # 3 components
    p.plot(n_components=3, switch_z=True)
