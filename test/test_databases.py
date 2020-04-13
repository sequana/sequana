from sequana import databases
import os
import glob
import pytest
import pytest_timeout

import requests

def network_available():
    try:
        res = requests.request("GET", "https://github.com", timeout=5)
        return True
    except requests.ConnectTimeout:
        return False
    except err:
        raise err

def test_eutils():
    if network_available():
        from sequana.databases import EUtilsTools
        et = EUtilsTools()
        res = et.accession_to_info("K01711.1")
        assert res['K01711.1']['gi'] == '331784'
        assert res['K01711.1']['taxid'] == '11234'


def test_database_download():
    if network_available():

        d = databases.ENADownload()
        d.download_viroid()
        for this in glob.glob('Viroid/*'):
            os.remove(this)
        os.rmdir("Viroid")

        d.download_list()
        # cleanup
        for key in d._metadata.keys():
            res = d._metadata[key]
            filename = res[0]
            os.remove(filename)

        d.download_accession('AF439431')
        os.remove("Custom/ENA_AF439431.fa")
        os.rmdir("Custom")

