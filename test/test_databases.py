from sequana import databases
import os
import glob
import pytest
import pytest_timeout

import requests
from . import test_dir


def network_available():
    try:
        res = requests.request("GET", "https://github.com", timeout=5)
        return True
    except requests.ConnectTimeout:
        return False
    except err:
        raise err


@pytest.mark.xfail(reason="too slow", method="thread")
@pytest.mark.timeout(10)
def test_eutils():
    if network_available():
        from sequana.databases import EUtilsTools

        et = EUtilsTools()
        res = et.accession_to_info("K01711.1")
        assert res["K01711.1"]["gi"] == 331784
        assert res["K01711.1"]["taxid"] == 11234
        sequence = et.get_fasta("K01711")


def test_NCBITaxonReader():
    # data test files were created using the taxons 11234 and 2697049

    # constructor with existing files
    n = databases.NCBITaxonReader(
        f"{test_dir}/data/names_filtered.dmp", f"{test_dir}/data/nodes_filtered.dmp"
    )

    # constructor zith no files
    n = databases.NCBITaxonReader()
    # required an init afterwards
    n.init(f"{test_dir}/data/names_filtered.dmp", f"{test_dir}/data/nodes_filtered.dmp")
    assert len(n.df_names) == 63
    assert len(n.df_nodes) == 22

    n.get_family(11234)

    # tempfile
    with open("output.dmp", "w") as fout:
        n.filter_names_dmp_file(output=fout.name, taxons=[11234])
        n.filter_nodes_dmp_file(output=fout.name, taxons=[11234])

    assert all(n.search("measles").taxon.values == [11234])

    n.get_average_name_per_taxon()
    assert n.get_number_taxon() == 22
    assert n.get_scientific_name(11234) == "Measles morbillivirus"
    assert n.get_taxon_from_scientific_name("Measles morbillivirus").values[0] == 11234


@pytest.mark.timeout(20)
@pytest.mark.xfail(reason="too slow")
def test_NCBIDownload(tmpdir):
    path = tmpdir.mkdir("temp").join("summary.txt")
    n = databases.NCBIDownload()
    try:
        filenames = []
        filenames = n.download_ncbi_refseq_release("mitochondrion")
    except Exception:
        assert False
    finally:
        for ff in filenames:
            os.remove(ff)
    n.download_assembly_report("fungi", output=path)



@pytest.mark.xfail(reason="too slow", method="thread")
@pytest.mark.timeout(10)
def test_ENADownload(tmpdir):
    path = tmpdir.mkdir("temp")
    n = databases.ENADownload()
    n.download_fasta("K01711.1", output_dir=str(path))
    n.download_fasta(["K01711.1", "dummy"], output_dir=str(path))
    n.download_fasta(
        f"{test_dir}/data/test_eutils_list_accession_number.txt", output_dir=str(path)
    )
