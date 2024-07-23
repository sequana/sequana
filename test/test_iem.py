from sequana.iem import SampleSheet as IEM
import glob

from . import test_dir


def test_iem():


    try:
        iem = IEM('dummy')
        assert False
    except:
        assert True


    for this in [
        "iem/wrong/test_expdesign_wrong.csv",
        "iem/good/test_expdesign_miseq_illumina_1.csv",
        "iem/good/test_expdesign_miseq_illumina2.csv",
    ]:

        filename = f"{test_dir}/data/{this}"

        iem = IEM(filename)
        iem.settings
        iem.name
        iem.samples
        iem.index_adapters
        iem.header
        iem.to_fasta("TEST")
        try:
            iem.validate()
        except SystemExit:
            pass
        iem.checker()
        print(iem)


def test_quick_fix(tmpdir):

    fout = tmpdir.join("fix.csv")
    e = IEM(f"{test_dir}/data/iem/wrong/test_expdesign_miseq_illumina_semicommas.csv")
    e.quick_fix(fout)
    e = IEM(fout)
    e.validate()


def test_warning():

    for filename in glob.glob(f"{test_dir}/data/iem/warning/*csv"):
        print(filename)
        e = IEM(filename)
        e.validate()
        checks = e.checker()
        errors = [1 for check in checks if check["status"] == "Error"]
        assert len(errors) == 0
        warnings = [1 for check in checks if check["status"] == "Warning"]
        assert len(warnings) > 0

def test_wrong():

    for filename in glob.glob(f"{test_dir}/data/iem/wrong/*csv"):

        e = IEM(f"{filename}")

        try:
            e.validate()
            assert False
        except SystemExit:
            checks = e.checker()
            errors = [1 for check in checks if check["status"] == "Error"]
            if sum(errors) == 0:
                assert False, filename



def test_iem_samplesheets():
    for filename in glob.glob(f"{test_dir}/data/iem/good/*csv"):
        e = IEM(f"{filename}")
        e.validate()
        e.version
        e.instrument
        e.header
        checks = e.checker()
        for check in checks:
            if check["status"] == "Error":
                assert False
