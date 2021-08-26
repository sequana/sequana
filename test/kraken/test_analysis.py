from sequana import sequana_config_path, sequana_data
from sequana import (
    KrakenAnalysis,
    KrakenDB,
    KrakenDownload,
    KrakenPipeline,
    KrakenSequential,
    KrakenResults,
)
import os
import tempfile
import pytest

from . import test_dir

database = sequana_config_path + os.sep + "kraken2_dbs/toydb"


@pytest.fixture
def download():
    kd = KrakenDownload()
    kd.download("toydb")


def test_download(tmpdir):
    # save in specific path
    p = tmpdir.mkdir("kr")
    kd = KrakenDownload(output_dir=str(p))
    kd.download("toydb")


def test_krakenDB(download):

    assert KrakenDB("toydb").__repr__() == "toydb"
    assert KrakenDB("toydb").version == "kraken2"
    # test from a full path
    k = KrakenDB(sequana_config_path + os.sep + "kraken2_dbs/toydb")
    k = KrakenDB(k)
    try:
        KrakenDB("dummy")
        assert False
    except IOError:
        assert True


@pytest.mark.xfail(reason="too slow or service may be down")
def test_run_kraken_analysis(download):

    file1 = sequana_data("Hm2_GTGAAA_L005_R1_001.fastq.gz", "data")
    file2 = sequana_data("Hm2_GTGAAA_L005_R2_001.fastq.gz", "data")
    kt = KrakenAnalysis([file1, file2], database=database, threads=1)
    kt.run()


@pytest.mark.xfail(reason="too slow or service may be down")
def test_kraken_sequential(download):
    p = tempfile.TemporaryDirectory()

    file1 = sequana_data("Hm2_GTGAAA_L005_R1_001.fastq.gz", "data")
    file2 = sequana_data("Hm2_GTGAAA_L005_R2_001.fastq.gz", "data")

    # Test paired data with 2 kraken1 DB
    kt = KrakenSequential(
        [file1, file2], [database, database], output_directory=p.name, force=True
    )
    kt.run()

    kt = KrakenSequential(
        file1, [database, database], output_directory=p.name, force=True
    )
    kt.run()

    p.cleanup()


@pytest.mark.xfail(reason="too slow or service may be down")
def test_kraken_results(download):
    k = KrakenResults(f"{test_dir}/data/test_kraken.out")
    df = k.plot(kind="pie")
    k.boxplot_classified_vs_read_length()
    print(df)

    df = k.plot(kind="barh")

    df = k.get_taxonomy_db(11234)
    assert 11234 in df.index

    from easydev import TempFile

    with TempFile() as fout:
        k.kraken_to_csv(fout.name, "toydb")
        k.kraken_to_json(fout.name, "toydb")
        k.kraken_to_krona(fout.name)
        k.to_js(fout.name)
    df = k.plot2(kind="pie")

    assert len(k.get_taxonomy_db([])) == 0
    k.histo_classified_vs_read_length()


@pytest.mark.xfail(reason="too slow or service may be down")
def test_kraken_pipeline(download):

    file1 = sequana_data("Hm2_GTGAAA_L005_R1_001.fastq.gz", "data")
    file2 = sequana_data("Hm2_GTGAAA_L005_R2_001.fastq.gz", "data")
    kp = KrakenPipeline([file1, file2], database=database, threads=1)
    kp.run()
