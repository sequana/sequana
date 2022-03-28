from sequana import sequana_data


def test_sequana_data():

    f1 = sequana_data("Institut_Pasteur.png")
    f2 = sequana_data("Institut_Pasteur.png", "images")
    assert f1 == f2


def test_sequana_data_star():
    # all files in a specific directory (a list)
    f1 = sequana_data("*", "images")
    assert isinstance(f1, list)
    assert 'Institut_Pasteur.png' in f1

    # all files (return a dict)
    f1 = sequana_data("*")
    assert isinstance(f1, dict)


def test_sequana_data():

    try:
        sequana_data()
        assert False
    except ValueError:
        assert True


    try:
        sequana_data("dummy")
        assert False
    except Exception:
        assert True


    try:
        sequana_data("dummy", "dummy")
        assert False
    except Exception:
        assert True


