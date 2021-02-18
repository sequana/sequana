from sequana.kegg import KEGGHelper



def test_kegg():
    k = KEGGHelper()
    results = k.search("lepto")
    print(results)

    results = k.search("9606")
    assert len(results) == 2

    from easydev import TempFile
    with TempFile() as fout:
        k.build_csv(fout.name, Nmax=2)

