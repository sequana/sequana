from sequana.assembly import BUSCO
from . import test_dir


def test_busco(tmpdir):

    outname = tmpdir.join('test.png')
    filename = f"{test_dir}/data/test_assembly/test_busco_full_table.tsv"
    b = BUSCO(filename)
    print(b)
    b.pie_plot(filename=outname)
    b.scatter_plot(filename=outname)
    assert b.score > 90
    assert b.score <91

def test_busco_v5():
    filename = f"{test_dir}/data/test_assembly/test_busco_full_table_v5.tsv"
    b = BUSCO(filename)
    print(b)
    b.pie_plot()
    b.scatter_plot()
    assert b.score == 100


def test_busco_core_genome(tmpdir):

    outname = tmpdir.join('test.fa')
    b = BUSCO(f"{test_dir}/data/test_assembly/busco_tiny.tsv")
    b.save_core_genomes(f"{test_dir}/data/test_assembly/data.contigs.fasta", outname)

    from sequana import FastA
    f = FastA(outname)
    assert len(f) == 5

