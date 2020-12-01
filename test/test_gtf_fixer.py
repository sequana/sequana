from sequana.gtf import GTFFixer
from sequana import sequana_data
from easydev import TempFile

prog = "sequana_gtf_fixer"

def test_input():
    filename = sequana_data('test_gtf_fixer.gtf')

    with TempFile() as fout:
        gtf = GTFFixer(filename)
        gtf.fix(fout.name)




