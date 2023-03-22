from sequana.gtf_fixer import GTFFixer

prog = "sequana_gtf_fixer"
from . import test_dir
def test_input(tmp_path):
    filename = f"{test_dir}/data/gtf/test_gtf_fixer.gtf"


    fout = tmp_path / "test.gtf"
    gtf = GTFFixer(filename)
    gtf.fix(fout)




