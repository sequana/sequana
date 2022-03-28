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
