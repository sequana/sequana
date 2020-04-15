
from sequana.gff3 import GFF3
from sequana import sequana_data
from easydev import TempFile

def test_gff():
    gff = GFF3(sequana_data('test_small.gff3'))
    df = gff.get_df()
    assert "telomere" in gff.get_types()


def test_gff_rnadiff():
    gff = GFF3(sequana_data('saccer3_truncated.gff'))
    df = gff.get_df()
    with TempFile() as fout:
        gff.create_files_for_rnadiff(fout.name)

        import pandas as pd
        df1 = pd.read_csv("{}_gene_lengths.tsv".format(fout.name), sep='\t')
        assert df1.Length.sum() == 31755 


