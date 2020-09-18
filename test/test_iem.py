


from sequana.iem import IEM
from sequana import sequana_data
from easydev import TempFile


def test_iem():

    for this in ["test_expdesign_wrong.csv",
        "test_expdesign_miseq_illumina_1.csv",
        "test_expdesign_miseq_illumina2.csv"]:

        filename = sequana_data(this)

        iem = IEM(filename)
        iem._scanner()
        iem.settings
        iem.name
        iem.header
        iem.to_fasta("TEST")
        try:iem.validate()
        except:pass


    with TempFile() as fh:

        e = IEM(sequana_data("test_expdesign_miseq_illumina_semicommas.csv"), tryme=True)
        e.quick_fix(fh.name)
        e = IEM(fh.name)


