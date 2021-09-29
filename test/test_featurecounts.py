import os
import glob

import sequana.featurecounts as fc
from sequana import sequana_data

from . import test_dir
test_dir = f"{test_dir}/data"

RNASEQ_DIR_0 = f"{test_dir}/featurecounts/rnaseq_0"
RNASEQ_DIR_1 = f"{test_dir}/featurecounts/rnaseq_1"
RNASEQ_DIR_2 = f"{test_dir}/featurecounts/rnaseq_2"
RNASEQ_DIR_undef = f"{test_dir}/featurecounts/rnaseq_undef"
RNASEQ_DIR_noconsensus = f"{test_dir}/featurecounts/rnaseq_noconsensus"
RNASEQ_DIR = f"{test_dir}/featurecounts/featurecounts_ex1"
RNASEQ_NEW_1= f"{test_dir}/featurecounts/new_rnaseq_output"


def __test_file():
    RNASEQ_DIR_0 = f"{test_dir}/featurecounts/rnaseq_0"
    assert os.path.isdir(RNASEQ_DIR_0)
    assert os.path.isdir(RNASEQ_DIR_0+"/sample1")
    assert os.path.isdir(RNASEQ_DIR_0+"/sample1/feature_counts_0")
    print(glob.glob(RNASEQ_DIR_0+"/sample1/feature_counts_0/*"))

def test_featurecounts():

    print( fc.get_most_probable_strand_consensus(RNASEQ_DIR_0, tolerance=0.1))

    assert fc.get_most_probable_strand_consensus(RNASEQ_DIR_0, tolerance=0.1)[0] == 0
    assert fc.get_most_probable_strand_consensus(RNASEQ_DIR_1, tolerance=0.1)[0] == 1
    assert fc.get_most_probable_strand_consensus(RNASEQ_DIR_2, tolerance=0.1)[0] == 2
    assert fc.get_most_probable_strand_consensus(RNASEQ_NEW_1, tolerance=0.1)[0] == 1

    assert (
        fc.get_most_probable_strand_consensus(RNASEQ_DIR_undef, tolerance=0.1)[0] == -1
    )

    # FIXME. how to handle the case where same amount of e.g 0 and 2 or 0 and 1
    # etc ? 
    #try:
    #    fc.get_most_probable_strand_consensus(RNASEQ_DIR_noconsensus, tolerance=0.1)
    #    assert False
    #except IOError:
    #    assert True

def test_multi_feature_counts():
    ff = fc.MultiFeatureCount(RNASEQ_DIR_0, 0.15)
    ff.get_most_probable_strand_consensus()
    ff.plot_strandness()

def test_feature_counts():
    fc1 = fc.FeatureCount(RNASEQ_DIR + "/all_features.out")
    fc2 = fc.FeatureCount(glob.glob(RNASEQ_DIR + "/*feature.out"))

    # we sort index to avoid error in python 3.6 travis
    assert all(fc1.df[fc1.df.columns] == fc2.df[fc1.df.columns])
    assert all(fc1.rnadiff_df[fc1.rnadiff_df.columns] == fc2.rnadiff_df[fc1.rnadiff_df.columns])











