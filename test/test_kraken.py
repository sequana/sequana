from sequana.kraken import *
from sequana import sequana_data, sequana_config_path
import os
import tempfile
import pytest


skiptravis = pytest.mark.skipif("TRAVIS_PYTHON_VERSION" in os.environ,
    reason="On travis")


#@skiptravis
def run_kraken_taxon():

    def download():
        kd = KrakenDownload()
        kd.download('toydb')
    database = sequana_config_path + os.sep + "kraken_toydb"
    if os.path.exists(database) is False:
        download()

    file1 = sequana_data("Hm2_GTGAAA_L005_R1_001.fastq.gz", "data")
    file2 = sequana_data("Hm2_GTGAAA_L005_R2_001.fastq.gz", "data")

    kt = KrakenAnalysis([file1, file2], database=database)
    kt.run()

    kt = KrakenAnalysis(file2, database=database)
    kt.run()

    p = tempfile.TemporaryDirectory()

    kt = KrakenHierarchical([file1, file2], [database, database],
            output_directory=p.name, force=True)
    kt.run()

    kt = KrakenHierarchical(file1, [database, database],
output_directory=p.name, force=True)
    kt.run()

    p.cleanup()



@skiptravis
def test_kraken_results():
    test_file = sequana_data("test_kraken.out", "testing")
    k = KrakenResults(test_file)
    df = k.plot(kind='pie')
    print(df)

    df = k.plot(kind='barh')

    df = k.get_taxonomy_biokit(11234)
    assert 11234 in df.index

    df = k.get_taxonomy_biokit("11234")
    assert 11234 in df.index


    df = k.get_taxonomy_biokit([11234])
    assert 11234 in df.index

    df = k.get_taxonomy_biokit(["11234"])
    assert 11234 in df.index


def test_download():
    kd = KrakenDownload()
    kd.download('toydb')

