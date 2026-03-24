from easydev import TempFile

from sequana.gff3 import GFF3

from . import test_dir


def test_wrong_input():
    try:
        gff = GFF3(f"{test_dir}/data/missing")
        gff.df
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
    g.get_attributes("gene")


def test_get_duplicated_attributes_per_genetic_type():

    g = GFF3(f"{test_dir}/data/gff/lenny.gff")
    g.get_duplicated_attributes_per_genetic_type()


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


def test_save_gff_filtered():
    gff = GFF3(f"{test_dir}/data/ecoli_truncated.gff")
    gff.get_seqid2size()


def test_contig_names():
    gff = GFF3(f"{test_dir}/data/ecoli_truncated.gff")
    names = gff.contig_names
    assert isinstance(names, list)
    assert len(names) > 0


def test_search():
    gff = GFF3(f"{test_dir}/data/ecoli_truncated.gff")
    result = gff.search("gene")
    assert len(result) > 0
    result_empty = gff.search("ZZZNOMATCH999")
    assert len(result_empty) == 0


def test_get_simplify_dataframe():
    gff = GFF3(f"{test_dir}/data/ecoli_truncated.gff")
    df = gff.get_simplify_dataframe()
    assert "genetic_type" in df.columns
    assert "seqid" in df.columns
    assert len(df) > 0


def test_add_CDS_and_mRNA(tmpdir):
    # Build a minimal GFF with gene features that add_CDS_and_mRNA can process
    gff_content = (
        "##gff-version 3\n"
        "chr1\ttest\tgene\t100\t200\t.\t+\t.\tgene_id=gene1;Name=gene1\n"
        "chr1\ttest\tgene\t300\t400\t.\t-\t.\tgene_id=gene2;Name=gene2\n"
    )
    infile = tmpdir.join("input.gff")
    infile.write(gff_content)
    outfile = tmpdir.join("output.gff")

    gff = GFF3(str(infile))
    gff.add_CDS_and_mRNA(str(outfile))

    content = outfile.read()
    assert "mRNA" in content
    assert "CDS" in content
    # Original gene lines preserved
    assert "gene1" in content
    assert "gene2" in content
