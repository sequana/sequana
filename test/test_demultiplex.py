from sequana.demultiplex import StatsFile
from sequana import sequana_data
from easydev import TempFile
import os

def test_stats_file():


    data = sequana_data("test_demultiplex_Stats.json")
    s = StatsFile(data)
    with TempFile() as fout:
        s.to_summary_reads(fout.name)
    with TempFile() as fout:
        s.barplot_summary(fout.name)
    with TempFile() as fout:
        s.barplot()
        for lane in s.get_data_reads().lane.unique():
             os.remove("lane{}_status.png".format(lane)) 

