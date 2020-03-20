
from sequana.gff3 import GFF3
from sequana import sequana_data


def test_gff():
    gff = GFF3(sequana_data('test_small.gff3'))
    df = gff.get_df()
    assert "telomere" in gff.get_types()
