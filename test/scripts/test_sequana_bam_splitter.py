from sequana.scripts import bam_splitter
from sequana import sequana_data
from tempfile import TemporaryDirectory
import os
import pytest

skiptravis = pytest.mark.skipif("TRAVIS_PYTHON_VERSION" in os.environ, reason="On travis")


prog = "sequana_bam_splitter"


def test_input():
    filename = sequana_data('test.bam')

    with TemporaryDirectory() as tmpdir:
        bam_splitter.main([prog, '--input', filename,
            "--output-directory" , tmpdir])

    with TemporaryDirectory() as tmpdir:
        bam_splitter.main([prog, '--input', filename,
            "--output-directory" , tmpdir, "--keep-unmapped", "--prefix", "test"])

    try: bam_splitter.main([prog, '--help'])
    except:pass
    try: bam_splitter.main([prog])
    except:pass
    try: bam_splitter.main()
    except:pass
    try: bam_splitter.main([prog, '--version'])
    except:pass


    #bam_splitter.main([prog, "--version"])

def test_output():
    with TemporaryDirectory() as tmpdir:

        prefix = tmpdir + "/test"
        M, U, F = bam_splitter._main(sequana_data("test.bam"), prefix,
            keep_unmapped=True)

        M, U, F = bam_splitter._main(sequana_data("test.bam"), prefix)
        from collections import Counter
        assert M == 934
        assert U == 66
        c = Counter(F)
        assert len(F) == 1000
        # ideally we should test all different flags. Here we test only a few of
        # them
        assert F.count(81) == 217
        assert F.count(73) == 1
        assert F.count(145) == 242
        assert F.count(163) == 229
        assert F.count(99) == 220

def test_sam_bam():
    with TemporaryDirectory() as tmpdir:
        prefix = tmpdir + "/test"
        M, U, F = bam_splitter._main(sequana_data("test_measles.bam"), prefix,
            keep_unmapped=True)

    with TemporaryDirectory() as tmpdir:
        prefix = tmpdir + "/test"
        M, U, F = bam_splitter._main(sequana_data("test_measles.sam"), prefix,
            keep_unmapped=True)

@pytest.mark.xfail(reason="cram flaky")
def test_cram():
    with TemporaryDirectory() as tmpdir:
        prefix = tmpdir + "/test"
        M, U, F = bam_splitter._main(sequana_data("test_measles.cram"), prefix,
            keep_unmapped=True)

def test_fasta():
    with TemporaryDirectory() as tmpdir:
        prefix = tmpdir + "/test"
        try:
            M, U, F = bam_splitter._main(sequana_data("test.fasta"), prefix,
                keep_unmapped=True)
        except ValueError:
            pass










