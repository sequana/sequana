from sequana import TRF

from . import test_dir

def test_trf():
    tt = TRF(f"{test_dir}/data/trf/test_trf1.dat")
    tt.hist_period_size()
    tt.hist_entropy()
    tt.hist_repet_by_sequence()
