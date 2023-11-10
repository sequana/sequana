from sequana import bedtools
from sequana.modules_report.coverage import CoverageModule
from sequana.utils import config


from .. import test_dir

def test_coverage_module(tmpdir):

    bed = bedtools.GenomeCov(f"{test_dir}/data/bed/JB409847.bed", 
        reference_file=f"{test_dir}/data/fasta/JB409847.fasta")
    c = bed[0]
    c.run(4001)

    directory = tmpdir.mkdir('test_coverage_module')
    config.output_dir = str(directory)
    config.sample_name = "JB409847"
    CoverageModule(bed)
