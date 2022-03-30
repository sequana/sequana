from sequana.cnv import CNVnator

from . import test_dir

def test_cnvnator():

    data = f"{test_dir}/data/txt/test_cnvnator.txt"
    c = CNVnator(data)
    c.plot("ENA|FN433596|FN433596.1")
    c.hist_event_size()
