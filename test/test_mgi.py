from sequana.mgi import MGI

from . import test_dir


def test_mgi():
    m = MGI(f"{test_dir}/data/mgi/test_mgi.fqStat.txt")
    m.plot_acgt()
    m.boxplot_quality()
