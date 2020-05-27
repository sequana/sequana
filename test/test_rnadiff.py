from sequana.rnadiff import RNADiffResults
from sequana import sequana_data




def test_rnadiff_onecond():
    RNADIFF_DIR = sequana_data("rnadiff") + "/rnadiff_onecond_1"
    r = RNADiffResults(RNADIFF_DIR)

    r.plot_count_per_sample()
    r.plot_percentage_null_read_counts()
    r.plot_volcano()
    r.pca()

def test_rnadiff_several_cond():
    RNADIFF_DIR = sequana_data("rnadiff") + "/rnadiff_analysis_1"
    r = RNADiffResults(RNADIFF_DIR)
    r.pca() 
