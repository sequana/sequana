from sequana.errors import SequanaException


def test_exception():
    try:
        raise SequanaException("test")
        assert False
    except:
        assert True
