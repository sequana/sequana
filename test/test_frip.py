from sequana.frip import FRiP


from . import test_dir

def test_frip():


    data = f"{test_dir}/data/FRiP/data.csv"
    fr = FRiP(data)
    fr.plot()
