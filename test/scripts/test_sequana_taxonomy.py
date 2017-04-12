from sequana.scripts import taxonomy
from sequana import sequana_data
import os
import pytest

prog = "sequana_taxonomy"


@pytest.fixture
def krakendb():
    # todo
    try:
        taxonomy.main([prog, '--download', 'toydb'])
    except TypeError: # Fails on travis so we download manually (appdirs returns
                      # none instead of the expected user config path 
        HOME = os.getenv('HOME')
        from sequana.misc import wget
        baseurl = "https://github.com/sequana/data/raw/master/kraken_toydb/"
        filenames = [
             "database.idx",
             "database.kdb",
             "taxonomy/names.dmp",
             "taxonomy/nodes.dmp"]
        for filename in filenames:
            from easydev import mkdirs
            mkdirs(HOME + os.sep + "database/taxonomy")
            wget(baseurl + os.sep + filename, 
                os.sep.join([HOME, "database", filename]))
    except SystemExit:
        pass


def test_analysis(krakendb):
    file1 = sequana_data("Hm2_GTGAAA_L005_R1_001.fastq.gz")
    file2 = sequana_data("Hm2_GTGAAA_L005_R2_001.fastq.gz")

    # Test that database must be provided
    try:
        df = taxonomy.main([prog, '--file1', file1])
        assert False
    except:
        assert True

    from tempfile import TemporaryDirectory
    directory = TemporaryDirectory()

    # If on travis and we could not load the database, use the local one
    # that must have been downloaded
    try:
        df = taxonomy.main([prog, '--file1', file1, "--database", "toydb",
            "--file2", file2, "--level", "INFO", "--output-directory", directory.name])
    except:
        HOME = os.getenv('HOME')
        database = os.sep.join([HOME, 'database'])
        df = taxonomy.main([prog, '--file1', file1, "--database", database,
            "--file2", file2,  "--output-directory", directory.name])
    from sequana import logger
    logger.info(directory.name)

def test_help():
    try:
        taxonomy.main([prog, '--help', '1>/tmp/out', '2>/tmp/err'])
        assert False
    except SystemExit:
        pass
    else:
        raise Exception

def test_wrong_db():
    try:
        df = taxonomy.main([prog, "--database", "dummy"])
        assert False
    except:
        assert True
