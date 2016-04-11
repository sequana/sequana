from sequana import snaketools, sequana_data
from sequana.snaketools import get_tagname, DOTParser
import os


testing = lambda x: sequana_data(x, "testing")

def test_dot_parser():
    s = DOTParser(testing("test_dag.dot"))
    s.add_urls()
    try:os.remove("test_dag.ann.dot")
    except:pass

def test_modules():
    assert "dag" in snaketools.modules.keys()
    assert snaketools.modules['dag'].endswith("Snakefile.dag")

def test_getcleanup_rules():
    filename =  snaketools.modules['fastq_sampling']
    try:
        snaketools.get_cleanup_rules(filename)
    except:
        pass

def test_snakemake_stats():
    # this is created using snakemake with the option "--stats stats.txt"
    s = snaketools.SnakeMakeStats(testing("test_snakemake_stats.txt"))
    s.plot()

def test_expanded_snakefile():
    filename = snaketools.modules['dag']
    try:
        s = snaketools.ExpandedSnakeFile(filename)
        s.expand()
    except ImportError:
        assert True
    #TODO remove the file just created (Snakefile.expanded)

def test_get_tagname():
    assert get_tagname("test") == "test"
    assert get_tagname("test.txt") == "test"
    assert get_tagname("prefix/test.txt") == "test"

def test_module():
    m = snaketools.Module('dag')
    m.description

def test_valid_config():
    s = snaketools.Module("fix_removal")
    config = snaketools.SequanaConfig(s.config)
