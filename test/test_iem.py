


from sequana.iem import IEM
from easydev import TempFile

from . import test_dir
def test_iem():

    for this in ["test_expdesign_wrong.csv",
        "test_expdesign_miseq_illumina_1.csv",
        "test_expdesign_miseq_illumina2.csv"]:

        filename = f"{test_dir}/data/csv/{this}"

        iem = IEM(filename)
        iem._scanner()
        iem.settings
        iem.name
        iem.header
        iem.to_fasta("TEST")
        try:iem.validate()
        except:pass
        print(iem)

    with TempFile() as fh:

        e = IEM(f"{test_dir}/data/csv/test_expdesign_miseq_illumina_semicommas.csv", tryme=True)
        e.quick_fix(fh.name)
        e = IEM(fh.name)


def test_iem_samplesheets():
    for filename in ['test_iem_samplesheet_hiseq_single.csv',
                     'test_iem_samplesheet_miseq_double.csv',
                     'test_iem_samplesheet_nextseq_single.csv',
                     'test_iem_samplesheet_iseq100.csv',
                     'test_iem_samplesheet_misnieq.csv']:
        iem = IEM(f"{test_dir}/data/csv/{filename}")
        iem.validate()
        iem.version
        iem.instrument
        iem.header

def test_old_samplesheet():
    for filename in ['test_expdesign_hiseq.csv', 'test_expdesign_hiseq_doubleindex.csv']:
        iem = IEM(f"{test_dir}/data/csv/{filename}")
        iem.validate()
        iem.version
        iem.instrument
