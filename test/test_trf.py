from sequana import TRF
from sequana import sequana_data



def test_trf():


    tt = TRF(sequana_data("test_trf1.dat"))
    tt.hist_period_size()
    tt.hist_entropy()
    tt.hist_repet_by_sequence()
