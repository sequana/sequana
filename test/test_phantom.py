from sequana.phantom import Phantom, PhantomPeaksReader

from . import test_dir


def test_PhantomPeakReader():

    p = PhantomPeaksReader(f"{test_dir}/data/phantom/spp.out")
    p.read(f"{test_dir}/data/phantom/spp.out")
    assert len(p.df) == 2
    p.plot_RSC()


def test_phantom():

    c = Phantom(binning=5)
    c.read_align(f"{test_dir}/data/phantom/data.align")
    results, df = c.run()
    stats = c.stats(results, df)
