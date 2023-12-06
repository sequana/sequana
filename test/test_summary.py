from sequana.summary import Summary

from . import test_dir


def test_summary():
    s = Summary("test2", sample_name="chr1", data={"mean": 1})
    assert s.data == {"mean": 1}
    assert s.version
    assert s.date
    d = s.as_dict()
    assert "name" in d
    assert "version" in d
    assert "data" in d
    assert "date" in d

    # test wrong constructor
    try:
        s = Summary("test")
        assert False
    except:
        assert True

    try:
        s = Summary("test", "test")
        assert False
    except:
        assert True

    # test data_description
    s = Summary("test2", data={"mean": 1})
    s.data_description = {"mean": "mean of the data set"}
    assert s.data_description == {"mean": "mean of the data set"}
    try:
        s.data_description = {"dummy": 1}
        assert False
    except:
        assert True

    from easydev import TempFile

    with TempFile(suffix=".json") as fh:
        s.to_json(fh.name)


def test_summary_from_json():

    s = Summary(f"{test_dir}/modules_report/data/test_summary_fastq_stats.json")
