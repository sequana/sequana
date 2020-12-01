from sequana.mgi import MGI
from sequana import sequana_data


def test_mgi():
    m = MGI(sequana_data("test_mgi.fqStat.txt"))
    m.plot_acgt()
    m.boxplot_quality()
