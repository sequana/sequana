from sequana import bedtools
from sequana.modules_report.coverage import CoverageModule, ChromosomeCoverageModule
from sequana.utils import config


from .. import test_dir

def test_coverage_module(tmpdir):

    directory = tmpdir.mkdir('test_coverage_module')
    config.output_dir = str(directory)
    config.sample_name = "JB409847"

    # html reports are created here for each chromosome
    bed = bedtools.SequanaCoverage(f"{test_dir}/data/bed/JB409847.bed", 
        reference_file=f"{test_dir}/data/fasta/JB409847.fasta")
    chrom = bed[0]
    results = chrom.run(4001)


    ROIs = results.get_rois()  
    ChromosomeCoverageModule(
        chrom,
        datatable=CoverageModule.init_roi_datatable(ROIs),
        options={"W": 4001, "k": 2, "ROIs": ROIs, "circular": False},
        command="test" , 
    )

    CoverageModule(bed)
