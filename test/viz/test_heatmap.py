from sequana.viz import Heatmap

from . import test_dir


def test_heatmap():
    filename = f"{test_dir}/data/test_heatmap.csv"
    import pandas as pd
    data = pd.read_csv(filename, skiprows=2, index_col=0)

    h = Heatmap(data)
    h.plot(cmap='hot')
    h.row_method= 'single'
    h.col_method= 'single'
    #    category_cols=[0,0,1,1], 
    #    category_rows=[0,1,2,0,0,1,2,2,2,1])

def test_doc_example():
    from sequana.viz import heatmap
    df = heatmap.get_heatmap_df()
    h = heatmap.Heatmap(df)
    h.category_column['A'] = 1
    h.category_column['B'] = 1
    h.category_column['C'] = 2
    h.category_column['D'] = 2

    h.category_row[2] = 2
    h.category_row[3] = 1
    h.category_row[0] = 1
    h.category_row[1] = 2

    h.plot()


def test_methods_and_metrics():
    from sequana.viz import heatmap
    df = heatmap.get_heatmap_df()
    
    h = heatmap.Heatmap(df, row_method="average", row_metric="jaccard",
        column_metric="jaccard", column_method="average")
    h.column_method = "average"
    h.column_metric = "jaccard"
    h.row_method = "average"
    h.row_metric = "jaccard"
    h.plot()

def test_misc():
    from sequana.viz import heatmap
    df = heatmap.get_heatmap_df()
    h = heatmap.Heatmap(df)
    h.plot(colorbar_position="top left")
    h.plot(colorbar_position="right")
    try:
        h.plot(colorbar_position="left")
        assert False
    except:
        assert True
    h.plot(gradient_span="min_to_max_centered")
    h.plot(gradient_span="only_max")
    h.plot(gradient_span="only_min")

def test_others():
    h = Heatmap(None, verbose=True)
    h = Heatmap(None, verbose=False)
    try:
        h = Heatmap(1, verbose=True)
        assert False
    except:
        assert True
