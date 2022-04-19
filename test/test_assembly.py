from sequana.assembly import BUSCO
from . import test_dir


def test_busco():
    filename = f"{test_dir}/data/tsv/test_busco_full_table.tsv"
    b = BUSCO(filename)
    print(b)
    b.pie_plot()
    b.scatter_plot()
    assert b.score > 90
    assert b.score <91

def test_busco_v5():
    filename = f"{test_dir}/data/tsv/test_busco_full_table_v5.tsv"
    b = BUSCO(filename)
    print(b)
    b.pie_plot()
    b.scatter_plot()
    assert b.score == 100
