import os

from easydev import TempFile

from sequana.demultiplex import StatsFile

from . import test_dir


def test_stats_file():

    data = f"{test_dir}/data/json/test_demultiplex_Stats.json"
    s = StatsFile(data)
    with TempFile() as fout:
        s.to_summary_reads(fout.name)
    with TempFile() as fout:
        s.barplot_summary(fout.name)
    with TempFile() as fout:
        s.barplot()
        for lane in s.get_data_reads().lane.unique():
            os.remove("lane{}_status.png".format(lane))

    data = f"{test_dir}/data/json/test_demultiplex_Stats_undetermined.json"
    s = StatsFile(data)
    with TempFile() as fout:
        s.to_summary_reads(fout.name)
    with TempFile() as fout:
        s.barplot_summary(fout.name)
    with TempFile() as fout:
        s.barplot()
        for lane in s.get_data_reads().lane.unique():
            os.remove("lane{}_status.png".format(lane))
