from sequana.scripts import fixR1R2
from sequana import sequana_data
from tempfile import TemporaryDirectory
import os


def test_fixR1R2():
    filename = sequana_data('Hm2_GTGAAA_L005_R2_001.fastq.gz')
    prog = "fixR1R2.py"
    with TemporaryDirectory() as tmpdir:
        import sys
        sys.argv = [prog, filename, tmpdir]
        fixR1R2.main()



    filename = sequana_data('Hm2_GTGAAA_L005_R1_001.fastq.gz')
    prog = "fixR1R2.py"
    with TemporaryDirectory() as tmpdir:
        import sys
        sys.argv = [prog, filename, tmpdir]
        fixR1R2.main()
        

    try:
        sys.argv = [prog]
        fixR1R2.main()
        assert False
    except:
        assert True
