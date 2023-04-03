from sequana.iem import IEM

from . import test_dir


def test_iem():

    for this in [
        "test_expdesign_wrong.csv",
        "test_expdesign_miseq_illumina_1.csv",
        "test_expdesign_miseq_illumina2.csv",
    ]:

        filename = f"{test_dir}/data/csv/{this}"

        iem = IEM(filename)
        iem._scanner()
        iem.settings
        iem.name
        iem.samples
        iem.index_adapters
        iem.header
        iem.to_fasta("TEST")
        try:
            iem.validate()
        except:
            pass
        print(iem)


def test_quick_fix(tmpdir):

    fout = tmpdir.join("fix.csv")

    e = IEM(f"{test_dir}/data/csv/test_expdesign_miseq_illumina_semicommas.csv", tryme=True)
    e.quick_fix(fout)
    e = IEM(fout)


def test_wrong():

    for filename in [
        "test_iem_samplesheet_duplicate.csv",
        "test_iem_samplesheet_duplicated_index1.csv",
        "test_iem_samplesheet_duplicated_index2.csv",
    ]:
        e = IEM(f"{test_dir}/data/csv/{filename}")
        try:
            e.validate()
            assert False
        except SystemExit:
            assert True


def test_iem_samplesheets():
    for filename in [
        "test_iem_samplesheet_hiseq_single.csv",
        "test_iem_samplesheet_miseq_double.csv",
        "test_iem_samplesheet_nextseq_single.csv",
        "test_iem_samplesheet_iseq100.csv",
        "test_iem_samplesheet_miniseq.csv",
    ]:
        iem = IEM(f"{test_dir}/data/csv/{filename}")
        iem.validate()
        iem.version
        iem.instrument
        iem.header
