from sequana.rnadiff import RNADiffResults, plot_venn
from sequana import sequana_data


RNADIFF_DIR = sequana_data("rnadiff") + "/rnadiff_analysis_1"


def test_rnadiff():
    RNADiffResults(RNADIFF_DIR)


def test_plot_venn():
    plot_venn(compa_list=[{"a", "b"}, {"a"}], labels=["t1", "t2"])
