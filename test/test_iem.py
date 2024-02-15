from sequana.iem import IEM

from . import test_dir


def test_iem():

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
    for filename in ["iem/wrong/test_iem_samplesheet_one_sample_pe.csv"]:
        e = IEM(f"{test_dir}/data/{filename}")
        e.validate()
        checks = e.checker()
        errors = [1 for check in checks if check["status"] == "Error"]
        assert len(errors) == 0


def test_wrong():

    for filename in [
        "iem/wrong/test_iem_samplesheet_duplicate.csv",
        "iem/wrong/test_iem_samplesheet_duplicated_index1.csv",
        "iem/wrong/test_iem_samplesheet_duplicated_index2.csv",
        "iem/wrong/test_no_data_section.csv",
        "iem/wrong/test_wrong_data_header.csv",
        "iem/wrong/test_wrong_data_header_index.csv",
        "iem/wrong/test_expdesign_wrong.csv",
        "iem/wrong/test_non_homogen_I7_length.csv",
        "iem/wrong/test_empty_data_section.csv",
        "iem/wrong/test_empty_data_section2.csv",
        "iem/wrong/test_iem_samplesheet_one_sample_se.csv",
    ]:
        e = IEM(f"{test_dir}/data/{filename}")
        try:
            e.validate()
            assert False
        except SystemExit:
            checks = e.checker()
            errors = [1 for check in checks if check["status"] == "Error"]
            if sum(errors) == 0:
                assert False, filename


def test_iem_samplesheets():
    for filename in [
        "iem/good/test_iem_samplesheet_hiseq_single.csv",
        "iem/good/test_iem_samplesheet_miseq_double.csv",
        "iem/good/test_iem_samplesheet_nextseq_single.csv",
        "iem/good/test_iem_samplesheet_iseq100.csv",
        "iem/good/test_iem_samplesheet_miniseq.csv",
        "iem/good/test_iem_samplesheet_one_sample_se.csv",
        "iem/good/test_iem_samplesheet_one_sample_pe.csv",
    ]:
        iem = IEM(f"{test_dir}/data/{filename}")
        iem.validate()
        iem.version
        iem.instrument
        iem.header
        checks = iem.checker()
        for check in checks:
            if check["status"] == "Error":
                assert False
