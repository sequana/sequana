def test_dendogram():
    from sequana.viz import dendogram, heatmap

    df = heatmap.get_heatmap_df()
    h = dendogram.Dendogram(df)
    h.plot()
