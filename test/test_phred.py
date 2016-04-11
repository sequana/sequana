from sequana import Quality
from sequana import phred


def test_quality():

    q = Quality('ABC')
    q.plot()
    assert q.mean_quality == 33

    q = phred.QualitySanger('ABC')
    q = phred.QualitySolexa('ABC')

def test_ascii_to_quality():
    assert phred.ascii_to_quality("!") == 0
    assert phred.ascii_to_quality("~") == 93


def test_quality_to_ascii():
    assert phred.quality_to_ascii(65) == "b"
    assert phred.quality_to_ascii(32) == "A"
    assert phred.quality_to_ascii(32, phred=64) == "`"


def test_quality_to_proba():
    assert phred.quality_to_proba_sanger(0) == 1
    assert phred.quality_to_proba_sanger(40) == 0.0001


def test_others():
    #sanger proba quality
    assert phred.proba_to_quality_sanger(0) == 93
    assert phred.proba_to_quality_sanger(0.0001) == 40
    assert phred.proba_to_quality_sanger(1) == 0
    assert phred.proba_to_quality_sanger(2) == 0


    # solexa proba quality
    assert phred.proba_to_quality_solexa(0) == 62
    assert abs(phred.proba_to_quality_solexa(0.0001) - 40) < 1e-3
    assert phred.proba_to_quality_solexa(0.99) == -5
    assert phred.proba_to_quality_solexa(2) == -5


    # solexa and sanger quality are similar. In this exampl, sanger ones
    # is slighlty larger
    assert phred.quality_solexa_to_quality_sanger(64) > 64.

    #
    # inverse here
    assert phred.quality_sanger_to_quality_solexa(64) < 64
    assert phred.quality_sanger_to_quality_solexa(64) > 63.99
