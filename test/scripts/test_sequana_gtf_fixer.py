from sequana.scripts import gtf_fixer
from sequana import sequana_data
from easydev import TempFile

prog = "sequana_gtf_fixer"

def test_input():
    filename = sequana_data('test_gtf_fixer.gtf')

    with TempFile() as fout:
        df = gtf_fixer.main([prog, '--input', filename, "--output", fout.name])




