from easydev import TempFile

from sequana.gff3 import GFF3

from . import test_dir


def test_wrong_input():
    try:
        gff = GFF3(f"{test_dir}/data/missing")
        assert False
    except IOError:
        assert True


def test_various_gff():
    gff = GFF3(f"{test_dir}/data/test_small.gff3")
    df = gff.df
    assert "telomere" in gff.features

    gff = GFF3(f"{test_dir}/data/ecoli_truncated.gff")
    assert gff.df.iloc[0].ID == "id0"

    gff = GFF3(f"{test_dir}/data/mm10_truncated.gff")
    assert gff.df.iloc[0].ID == "chromosome:1"

    gff = GFF3(f"{test_dir}/data/hg38_truncated_gtf.gff")
    assert gff.df.iloc[0].gene_id == "ENSG00000223972"

    assert gff.clean_gff_line_special_characters("A%09A") == "A\tA"

    gff = GFF3(f"{test_dir}/data/gff/lenny.gff")
    df = gff.df

    gff = GFF3(f"{test_dir}/data/gff/Ld1S.gff")
    df = gff.df


def test_process_attributes():
    gff = GFF3(f"{test_dir}/data/mm10_truncated.gff")
    res = gff._process_attributes("ID=YAL058W;Name=YAL058W")
    assert len(res) == 2
    res = gff._process_attributes("ID YAL058W;Name YAL058W")
    assert len(res) == 2


def test_transcript_to_gene():
    gff = GFF3(f"{test_dir}/data/ecoli_truncated.gff")
    gff.transcript_to_gene_mapping(attribute="Name")


def test_read_and_save_selected_features(tmpdir):
    tmpfile = tmpdir.join("test.gff")
    gff = GFF3(f"{test_dir}/data/ecoli_truncated.gff")
    gff.read_and_save_selected_features(tmpfile)


def test_get_feature_dict():
    gff = GFF3(f"{test_dir}/data/ecoli_truncated.gff")
    gff.features
    gff.get_features_dict()


def test_attributes(tmpdir):
    g = GFF3(f"{test_dir}/data/gff/lenny.gff")
    assert g.attributes


def test_to_bed(tmpdir):
    outname = tmpdir.join("test.bed")
    gff = GFF3(f"{test_dir}/data/ecoli_truncated.gff")
    gff.to_bed(outname, "gene")


def test_to_fasta(tmpdir):
    outname = tmpdir.join("test.fasta")
    gff = GFF3(f"{test_dir}/data/gff/ecoli_MG1655.gff")
    gff.to_fasta(f"{test_dir}/data/fasta/ecoli_MG1655.fa", outname)


def test_gff_to_gtf():
    gff = GFF3(f"{test_dir}/data/saccer3_truncated.gff")
    with TempFile() as fout:
        df = gff.to_gtf(fout.name)


def test_save_gff_filtered():
    gff = GFF3(f"{test_dir}/data/saccer3_truncated.gff")
    with TempFile() as fout:
        df = gff.save_gff_filtered(filename=fout.name)
