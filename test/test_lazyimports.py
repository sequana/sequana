


def test_lazy(mocker):
    import sys
    import sequana.lazy as lazy
    import sequana.lazyimports as lazyimports
    import importlib as imp

    li = lazyimports.LazyImport('os')
    li
    # we import sphinx now, and reload the module so enter in the case where
    # sphinx is loaded. No need to install sphinx, we can mock it
    sys.modules['sphinx'] = "test"
    imp.reload(lazyimports)
    try:
        assert lazy.enabled() == False
        li = lazyimports.LazyImport("os")
        li.path
    except Exception as err:
        raise(err)
    finally:
        pass
    
