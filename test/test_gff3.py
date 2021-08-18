
from sequana.gff3 import GFF3
from easydev import TempFile

from . import test_dir

def test_gff():
    gff = GFF3(f"{test_dir}/data/test_small.gff3")
    df = gff.df
    assert "telomere" in gff.features


    gff = GFF3(f'{test_dir}/data/ecoli_truncated.gff')
    assert gff.df.iloc[0].ID == "id0"

    gff = GFF3(f'{test_dir}/data/mm10_truncated.gff')
    assert gff.df.iloc[0].ID == "chromosome:1"

    gff = GFF3(f'{test_dir}/data/hg38_truncated_gtf.gff')
    assert gff.df.iloc[0].gene_id == "ENSG00000223972"

def test_process_attributes():
    gff = GFF3(f'{test_dir}/data/mm10_truncated.gff')
    res = gff._process_attributes("ID=YAL058W;Name=YAL058W")
    assert len(res) == 2
    res = gff._process_attributes("ID YAL058W;Name YAL058W")
    assert len(res) == 2


def test_gff_rnadiff():
    gff = GFF3(f'{test_dir}/data/saccer3_truncated.gff')
    df = gff.df
    gff.get_duplicated_attributes_per_type()
    with TempFile() as fout:
        gff.create_files_for_rnadiff(fout.name)

        import pandas as pd
        df1 = pd.read_csv("{}_gene_lengths.tsv".format(fout.name), sep='\t')
        # changed in 21/07/2021 from 31755 to 29199.
        # not clear why the sum was 31755 before. 291999 check manually in the
        # GFF file
        assert df1.Length.sum() == 29199


    with TempFile() as fout:
        gff.save_annotation_to_csv(fout.name)

def test_gff_to_gtf():
    gff = GFF3(f'{test_dir}/data/saccer3_truncated.gff')
    with TempFile() as fout:
        df = gff.to_gtf(fout.name)

def test_save_gff_filtered():
    gff = GFF3(f'{test_dir}/data/saccer3_truncated.gff')
    with TempFile() as fout:
        df = gff.save_gff_filtered(filename=fout.name)
