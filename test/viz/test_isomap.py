

def test_isomap():
    from sequana.viz.isomap import Isomap
    from sequana import sequana_data
    import pandas as pd

    data = sequana_data("test_pca.csv")
    df = pd.read_csv(data)
    df = df.set_index("Id")
    p = Isomap(df, colors={
        "A1": 'r', "A2": 'r', 'A3': 'r',
        "B1": 'b', "B2": 'b', 'B3': 'b'})
    p.plot(n_components=2, switch_y=True)
    p.plot(n_components=2, switch_x=True)
    p.plot(n_components=3, switch_z=True)

