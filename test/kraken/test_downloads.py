
import pytest

from sequana import KrakenDownload 
from . import test_dir


def test_download(tmpdir):
    # save in specific path
    p = tmpdir.mkdir("kr")
    kd = KrakenDownload(output_dir=str(p))
    kd.download("toydb")

    # redownload on purpose
    kd.download("toydb")

