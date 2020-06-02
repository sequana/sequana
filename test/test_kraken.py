from sequana.kraken import *
from sequana import sequana_data, sequana_config_path
import os
import tempfile
import pytest


skiptravis = pytest.mark.skipif("TRAVIS_PYTHON_VERSION" in os.environ,
    reason="On travis")

try:
    kd = KrakenDownload()
    kd.download('toydb')
    kd.download('toydb2')
except:
    pass

@pytest.mark.xfail
def test_run_kraken_taxon():

    database = sequana_config_path + os.sep + "kraken_toydb"
    if os.path.exists(database) is False:
        download()

    file1 = sequana_data("Hm2_GTGAAA_L005_R1_001.fastq.gz", "data")
    file2 = sequana_data("Hm2_GTGAAA_L005_R2_001.fastq.gz", "data")
    kt = KrakenAnalysis([file1, file2], database=database, threads=1)
    kt.run()


@pytest.mark.xfail
def test_kraken_sequential():
    database = sequana_config_path + os.sep + "kraken_toydb"
    database2 = sequana_config_path + os.sep + "kraken2_dbs/toydb_nomasking"
    p = tempfile.TemporaryDirectory()

    file1 = sequana_data("Hm2_GTGAAA_L005_R1_001.fastq.gz", "data")
    file2 = sequana_data("Hm2_GTGAAA_L005_R2_001.fastq.gz", "data")

    # Test paired data with 2 kraken1 DB
    kt = KrakenSequential([file1, file2], [database, database2],
            output_directory=p.name, force=True)
    kt.run()

    kt = KrakenSequential(file1, [database, database],
        output_directory=p.name, force=True)
    kt.run()

    p.cleanup()


@pytest.mark.xfail
def test_kraken_results():
    test_file = sequana_data("test_kraken.out", "testing")
    k = KrakenResults(test_file)
    df = k.plot(kind='pie')
    k.boxplot_classified_vs_read_length()
    print(df)

    df = k.plot(kind='barh')

    df = k.get_taxonomy_db(11234)
    assert 11234 in df.index

    from easydev import TempFile
    with TempFile() as fout:
        k.kraken_to_csv(fout.name, "toydb")
        k.kraken_to_json(fout.name, "toydb")
        k.kraken_to_krona(fout.name )
        k.to_js(fout.name)
    df = k.plot2(kind='pie')
    

@pytest.mark.xfail
def test_kraken_pipeline():
    
    from sequana import KrakenPipeline
    file1 = sequana_data("Hm2_GTGAAA_L005_R1_001.fastq.gz", "data")
    file2 = sequana_data("Hm2_GTGAAA_L005_R2_001.fastq.gz", "data")
    database = sequana_config_path + os.sep + "kraken_toydb"
    kp = KrakenPipeline([file1, file2], database=database, threads=1)
    kp.run()



def test_mkr():
    from sequana.kraken import MultiKrakenResults
    mkr = MultiKrakenResults([sequana_data("test_kraken_multiple_1.csv"),
            sequana_data('test_kraken_multiple_1.csv')])
    mkr.plot_stacked_hist(kind="bar")           
    mkr.plot_stacked_hist(kind="barh")           


def test_mkr2():
    from sequana.kraken import MultiKrakenResults2
    mkr = MultiKrakenResults2([sequana_data("test_kraken_mkr2_summary_1.json"),
            sequana_data('test_kraken_mkr2_summary_2.json')])
    mkr.plot_stacked_hist()
