from sequana.cnv import CNVnator
from sequana import sequana_data


def test_cnvnator():

    data = sequana_data("test_cnvnator.txt")
    c = CNVnator(data)
    c.plot("ENA|FN433596|FN433596.1")
    c.hist_event_size()
