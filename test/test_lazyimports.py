# pytest Works locally but not at the package top level...


# This works but intervene with test_sphinxext somehow

def test_lazy(mocker):
    import sys
    import sequana.lazy as lazy
    import sequana.lazyimports as lazyimports
    import importlib as imp

    li = lazyimports.LazyImport('os')
    li
    # we import sphinx now, and reload the module so enter in the case where
    # sphinx is loaded
    import sphinx
    imp.reload(lazyimports)
    try:
        assert lazy.enabled() == False
        li = lazyimports.LazyImport("os")
        li.path
    except Exception as err:
        raise(err)
    finally:
        #del sys.modules['sphinx']
        #imp.reload(lazyimports)
        #assert lazy.enabled() == True
        pass
