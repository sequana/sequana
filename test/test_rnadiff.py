from sequana.rnadiff import RNADiffResults
from sequana import sequana_data




def test_rnadiff_onefolder():
    RNADIFF_DIR = sequana_data("rnadiff") + "/rnadiff_onecond_1"
    r = RNADiffResults(RNADIFF_DIR)
    r.plot_count_per_sample()
    r.plot_percentage_null_read_counts()
    r.plot_volcano()
    r.plot_pca()
    r.plot_mds()
    r.plot_isomap()

    try:
        r.get_cond_from_sample("fake")
        assert False
    except:
        assert True

def test_rnadiff_onefile():
    RNADIFF_DIR = sequana_data("rnadiff") + "/rnadiff_onecond_1"
    r = RNADiffResults(RNADIFF_DIR + "/tables/B3789-v1.surexpvsref.complete.xls")
    r.plot_count_per_sample()
    r.summary()

def test_rnadiff_several_cond():
    RNADIFF_DIR = sequana_data("rnadiff") + "rnadiff_fake_two_conditions"
    try:
        r = RNADiffResults(RNADIFF_DIR)
        assert False
    except TypeError:
        assert True
