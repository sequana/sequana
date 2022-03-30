from sequana.modules_report.fastq_stats import FastQStatsModule
import shutil
from sequana.utils import config



from . import test_dir
filename = f"{test_dir}/data/test_summary_fastq_stats.json"

def test_report(tmpdir):
    directory = tmpdir.mkdir('test_module')
    shutil.copy(filename, str(directory))
    config.output_dir = str(directory) 

    #
    #report = FastQStatsModule(tmpdir, "dummy_link_fastqc")

    #report.create_report()
