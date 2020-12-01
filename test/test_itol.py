from easydev import TempFile
from sequana import ITOL, sequana_data




# FIXME need a valid IP and API
def _test_itol_format():


    itol = ITOL(sequana_data("test_itol_basic.tree.txt"))
    itol.upload()
    for frmt in [".png", ".pdf", ".svg", ".eps"]:
        with TempFile(suffix=frmt) as fout:
            itol.export(fout.name)

    try:
        itol.export("dummy.gg")
        assert False
    except:
        assert True
    


def test_itol_status():
    # an error
    itol = ITOL(sequana_data("test_itol_basic.tree.txt"))

    with TempFile(suffix="png") as fout:
        # let us forget the upload    itol.upload()
        try:
            itol.export(fout.name)
            assert False
        except:
            assert True


