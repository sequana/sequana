from sequana.gtf import GTFFixer
from easydev import TempFile

prog = "sequana_gtf_fixer"
from . import test_dir
def test_input():
    filename = f"{test_dir}/data/gtf/test_gtf_fixer.gtf"

    with TempFile() as fout:
        gtf = GTFFixer(filename)
        gtf.fix(fout.name)




