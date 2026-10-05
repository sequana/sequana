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
        gff.save_gff_filtered(filename=fout.name)


def test_get_seqid2size_ecoli_truncated():
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


def test_get_intergenic_regions():
    gff = GFF3(f"{test_dir}/data/gff/ecoli_MG1655.gff")
    regions = gff.get_intergenic_regions()
    assert "seqid" in regions.columns
    assert "start" in regions.columns
    assert "stop" in regions.columns
    assert len(regions) > 0
    # every intergenic region must start before it ends
    assert (regions["stop"] >= regions["start"]).all()


def test_save_annotation_to_csv(tmpdir):
    gff = GFF3(f"{test_dir}/data/gff/ecoli_MG1655.gff")
    outfile = str(tmpdir.join("annotations.csv"))
    gff.save_annotation_to_csv(outfile)
    with open(outfile) as f:
        content = f.read()
    assert "seqid" in content
    assert "genetic_type" in content


def test_read_and_save_selected_features_only_requested_types(tmpdir):
    gff = GFF3(f"{test_dir}/data/gff/ecoli_MG1655.gff")
    outfile = str(tmpdir.join("genes_only.gff"))
    gff.read_and_save_selected_features(outfile, features=["gene"])
    with open(outfile) as f:
        lines = [line for line in f if line.strip()]
    assert len(lines) == 14
    assert all("\tgene\t" in line for line in lines)


def test_get_duplicated_attributes_per_genetic_type_returns_gene_column():
    gff = GFF3(f"{test_dir}/data/gff/ecoli_MG1655.gff")
    result = gff.get_duplicated_attributes_per_genetic_type()
    assert "gene" in result.columns

    result2 = gff.get_duplicated_attributes_per_genetic_type2()
    assert "gene" in result2.columns


def test_get_seqid2size():
    gff = GFF3(f"{test_dir}/data/gff/ecoli_MG1655.gff")
    sizes = gff.get_seqid2size()
    assert sizes["NC_000913.3"] == 4641652


def test_transcript_to_gene_mapping_error_on_missing_attribute():
    gff = GFF3(f"{test_dir}/data/gff/ecoli_MG1655.gff")
    # this GFF has no transcript_id attribute -> KeyError from the missing column
    try:
        gff.transcript_to_gene_mapping()
        assert False, "expected KeyError"
    except KeyError:
        pass


def test_process_attributes_missing_separator_raises_bad_file_format():
    from sequana.errors import BadFileFormat

    gff = GFF3(f"{test_dir}/data/ecoli_truncated.gff")
    try:
        gff._process_attributes("not_a_valid_attribute_string_without_separator")
        assert False, "expected BadFileFormat"
    except BadFileFormat:
        pass


def test_to_pep_not_implemented():
    gff = GFF3(f"{test_dir}/data/ecoli_truncated.gff")
    try:
        gff.to_pep("dummy.fasta", "dummy_out.fasta")
        assert False, "expected NotImplementedError"
    except NotImplementedError:
        pass


def test_add_directon_index_and_cluster_names_to_bed_not_implemented():
    gff = GFF3(f"{test_dir}/data/gff/ecoli_MG1655.gff")
    try:
        gff.add_directon_index()
        assert False, "expected NotImplementedError"
    except NotImplementedError:
        pass

    try:
        gff.cluster_names_to_bed()
        assert False, "expected NotImplementedError"
    except NotImplementedError:
        pass


def test_add_intergenic_regions():
    gff = GFF3(f"{test_dir}/data/gff/ecoli_MG1655.gff")
    n_before = len(gff.df)
    gff.add_intergenic_regions()
    assert len(gff.df) > n_before
    assert "region" in set(gff.df["genetic_type"])

    # second call is a documented no-op: the dataframe must not grow again
    n_after_first = len(gff.df)
    gff.add_intergenic_regions()
    assert len(gff.df) == n_after_first


def test_save_as_gff(tmpdir):
    gff = GFF3(f"{test_dir}/data/gff/ecoli_MG1655.gff")
    outfile = str(tmpdir.join("saved.gff"))
    gff.save_as_gff(outfile)
    with open(outfile) as f:
        content = f.read()
    assert "gene-b0001" in content
    assert "sorting seqid" in content


def test_save_as_gff_after_add_intergenic_regions(tmpdir):
    gff = GFF3(f"{test_dir}/data/gff/ecoli_MG1655.gff")
    gff.add_intergenic_regions()
    outfile = str(tmpdir.join("saved.gff"))
    gff.save_as_gff(outfile)
    with open(outfile) as f:
        content = f.read()
    assert "added intergenic region" in content


def test_save_gff_filtered_only_genes(tmpdir):
    gff = GFF3(f"{test_dir}/data/gff/ecoli_MG1655.gff")
    outfile = str(tmpdir.join("filtered.gff"))
    gff.save_gff_filtered(outfile, features=["gene"])
    with open(outfile) as f:
        content = f.read()
    assert "gene-b0001" in content
    assert "cds-NP_414542.1" not in content


def test_add_regions_and_save_gff(tmpdir):
    gff = GFF3(f"{test_dir}/data/gff/ecoli_MG1655.gff")
    outfile = str(tmpdir.join("with_regions.gff"))
    gff.add_regions_and_save_gff(outfile)
    with open(outfile) as f:
        content = f.read()
    assert "NC_000913.3" in content


def test_get_PTU_and_directons():
    """get_PTU() exercises _filter_coding_genes(), _compute_directons(), and
    the cached .directons property in one pass on a real bacterial GFF."""
    gff = GFF3(f"{test_dir}/data/gff/ecoli_MG1655.gff")
    ptu = gff.get_PTU()
    assert set(["chromosome", "start", "stop", "strand", "length"]) <= set(ptu.columns)
    assert len(ptu) > 0

    # directons property is cached: second access must return the same object
    first = gff.directons
    second = gff.directons
    assert first is second


def test_get_ssr():
    gff = GFF3(f"{test_dir}/data/gff/ecoli_MG1655.gff")
    ssr = gff._get_ssr()
    assert set(["type", "chromosome", "start", "stop"]) <= set(ssr.columns)
    assert set(ssr["type"]).issubset({"dSSR", "cSSR", "other"})


def test_directon_to_bed(tmpdir):
    gff = GFF3(f"{test_dir}/data/gff/ecoli_MG1655.gff")
    outfile = str(tmpdir.join("directons.bed"))
    gff.directon_to_bed(outfile)
    with open(outfile) as f:
        content = f.read()
    assert "NC_000913.3" in content
    # BED9: strand colours present
    assert "255,0,0" in content or "0,0,255" in content


def test_is_tRNA_or_rRNA():
    gff = GFF3(f"{test_dir}/data/gff/ecoli_MG1655.gff")
    assert gff._is_tRNA_or_rRNA("tRNA") is True
    assert gff._is_tRNA_or_rRNA("tRNA-Leu") is True
    assert gff._is_tRNA_or_rRNA("28S_rRNA") is True
    assert gff._is_tRNA_or_rRNA("thrA") is False
    assert gff._is_tRNA_or_rRNA(None) is False

    # deprecated inverted-semantics wrapper
    assert gff.is_tRNA_or_ribosomal("tRNA") is False
    assert gff.is_tRNA_or_ribosomal("thrA") is True
