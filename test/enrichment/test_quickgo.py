import pandas as pd

from sequana.enrichment.quickgo import QuickGOGraph

from . import test_dir


def test_get_graph():
    u = QuickGOGraph()
    IDs = ["GO:0016829", "GO:0051287", "GO:0000166", "GO:0022857", "GO:0051287"]
    u._get_graph(IDs, ["BP", "MF", "CC"])


def test_save_chart(tmpdir):
    outpng = tmpdir.join("test.png")
    u = QuickGOGraph()
    df = pd.read_csv(f"{test_dir}/data/test_quickgo_save_chart.csv")
    u.save_chart(df, outpng)


def test_get_go_description():

    u = QuickGOGraph()
    u.get_go_description(["GO:0016829", "GO:0051287", "GO:0000166"])
