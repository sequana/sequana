"""Tests for sequana.hmmtools (hmmscan domtblout parsing)."""
import pytest

from sequana.hmmtools import PfamDomtblout

# A minimal, well-formed hmmscan --domtblout fixture: 22 whitespace-separated
# fixed fields followed by a free-text description (which may itself contain
# spaces, hence the special parsing in PfamDomtblout.read()).
DOMTBLOUT_CONTENT = (
    "# hmmscan domtblout\n"
    "PF00001.21 accA 200 gene1.t1 accB 150 1.2e-30 100.5 0.0 1 1 1.2e-32 3.4e-30 99.5 0.0 2 190 10 155 8 158 0.95 7tm_1 receptor domain\n"
    "PF00002.15 accC 300 gene1.t1 accD 250 5.0e-10 50.0 0.0 1 1 5.0e-12 6.0e-10 45.0 0.0 5 280 20 200 15 210 0.80 second domain\n"
    "PF00003.10 accE 100 gene2.t1 accF 90 1.0e-05 20.0 0.0 1 1 1.0e-06 2.0e-05 15.0 0.0 1 90 5 85 3 88 0.70 third domain here\n"
)

AUGUSTUS_TRANSCRIPT_GFF = (
    "##gff-version 3\n"
    "chr1\tAUGUSTUS\tgene\t1000\t2000\t.\t+\t.\tID=gene1.t1\n"
    "chr1\tAUGUSTUS\tgene\t3000\t4000\t.\t+\t.\tID=gene2.t1\n"
)

AUGUSTUS_GENE_GFF = (
    "##gff-version 3\n"
    "chr1\tAUGUSTUS\tgene\t1000\t2000\t.\t+\t.\tID=gene1\n"
    "chr1\tAUGUSTUS\tgene\t3000\t4000\t.\t+\t.\tID=gene2\n"
    "chr1\tAUGUSTUS\tgene\t5000\t6000\t.\t+\t.\tID=gene3\n"
)


@pytest.fixture
def domtblout_file(tmpdir):
    path = tmpdir.join("pfam.domtblout")
    path.write(DOMTBLOUT_CONTENT)
    return str(path)


@pytest.fixture
def augustus_transcript_gff(tmpdir):
    path = tmpdir.join("augustus_transcript.gff")
    path.write(AUGUSTUS_TRANSCRIPT_GFF)
    return str(path)


@pytest.fixture
def augustus_gene_gff(tmpdir):
    path = tmpdir.join("augustus_gene.gff")
    path.write(AUGUSTUS_GENE_GFF)
    return str(path)


class TestPfamDomtbloutRead:
    def test_read_parses_all_rows(self, domtblout_file):
        p = PfamDomtblout(domtblout_file)
        p.read()
        assert len(p.df) == 3

    def test_read_skips_comment_lines(self, domtblout_file):
        p = PfamDomtblout(domtblout_file)
        p.read()
        assert not p.df["target_name"].str.startswith("#").any()

    def test_read_column_names(self, domtblout_file):
        p = PfamDomtblout(domtblout_file)
        p.read()
        assert p.df.columns[0] == "target_name"
        assert p.df.columns[-1] == "description"
        assert len(p.df.columns) == 23

    def test_read_handles_multiword_description(self, domtblout_file):
        p = PfamDomtblout(domtblout_file)
        p.read()
        assert p.df.iloc[0]["description"] == "7tm_1 receptor domain"
        assert p.df.iloc[2]["description"] == "third domain here"

    def test_read_fixed_field_values(self, domtblout_file):
        p = PfamDomtblout(domtblout_file)
        p.read()
        row0 = p.df.iloc[0]
        assert row0["target_name"] == "PF00001.21"
        assert row0["query_name"] == "gene1.t1"
        assert row0["i_evalue"] == "3.4e-30"


class TestPfamDomtbloutToGff:
    def test_to_gff_before_read_raises(self, tmpdir):
        p = PfamDomtblout("dummy.domtblout")
        with pytest.raises(ValueError):
            p.to_gff(str(tmpdir.join("out.gff")))

    def test_to_gff_best_hit_only(self, domtblout_file, tmpdir):
        p = PfamDomtblout(domtblout_file)
        p.read()
        outfile = str(tmpdir.join("out.gff"))
        p.to_gff(outfile, best_hit=True)

        with open(outfile) as f:
            lines = [line for line in f if not line.startswith("#")]
        # one line per query_name (gene1.t1 and gene2.t1) when best_hit=True
        assert len(lines) == 2

    def test_to_gff_all_hits(self, domtblout_file, tmpdir):
        p = PfamDomtblout(domtblout_file)
        p.read()
        outfile = str(tmpdir.join("out.gff"))
        p.to_gff(outfile, best_hit=False)

        with open(outfile) as f:
            lines = [line for line in f if not line.startswith("#")]
        # gene1.t1 has 2 hits, gene2.t1 has 1 hit -> 3 total
        assert len(lines) == 3

    def test_to_gff_with_augustus_translates_coordinates(self, domtblout_file, augustus_transcript_gff, tmpdir):
        p = PfamDomtblout(domtblout_file)
        p.read()
        outfile = str(tmpdir.join("out.gff"))
        p.to_gff(outfile, augustus_gff=augustus_transcript_gff, best_hit=True)

        with open(outfile) as f:
            lines = [line for line in f if not line.startswith("#") and line.strip()]
        assert len(lines) == 2
        fields = lines[0].split("\t")
        # seqid replaced by the augustus chromosome, start/stop offset by gene start (1000)
        assert fields[0] == "chr1"
        assert fields[3] == "1010"  # 1000 + ali_start(10)
        assert fields[4] == "1155"  # 1000 + ali_end(155)

    def test_to_gff_augustus_missing_match_is_skipped_not_raised(self, domtblout_file, tmpdir):
        """Regression test: previously this raised an unhandled IndexError
        when a query_name had no matching record in the Augustus GFF."""
        p = PfamDomtblout(domtblout_file)
        p.read()

        # augustus GFF that does NOT contain gene1.t1 / gene2.t1 at all
        empty_match_gff = tmpdir.join("no_match.gff")
        empty_match_gff.write("##gff-version 3\nchr1\tAUGUSTUS\tgene\t1\t100\t.\t+\t.\tID=unrelated_gene\n")

        outfile = str(tmpdir.join("out.gff"))
        p.to_gff(outfile, augustus_gff=str(empty_match_gff), best_hit=True)  # must not raise

        with open(outfile) as f:
            lines = [line for line in f if not line.startswith("#") and line.strip()]
        assert lines == []


class TestPfamDomtbloutAnnotateGff:
    def test_annotate_gff_adds_pfam_annotation_to_matching_genes(self, domtblout_file, augustus_gene_gff, tmpdir):
        p = PfamDomtblout(domtblout_file)
        p.read()
        outfile = str(tmpdir.join("annotated.gff"))
        p.annotate_gff(augustus_gene_gff, outfile)

        with open(outfile) as f:
            content = f.read()

        assert "ID=gene1;annot=7tm_1 receptor domain;pfamaccA" in content
        assert "ID=gene2;annot=third domain here;pfamaccE" in content

    def test_annotate_gff_uses_none_for_unmatched_genes(self, domtblout_file, augustus_gene_gff, tmpdir):
        p = PfamDomtblout(domtblout_file)
        p.read()
        outfile = str(tmpdir.join("annotated.gff"))
        p.annotate_gff(augustus_gene_gff, outfile)

        with open(outfile) as f:
            content = f.read()

        # gene3 has no hmmscan hit at all
        assert "ID=gene3;annot=none;pfamnone" in content

    def test_annotate_gff_preserves_non_gene_lines(self, domtblout_file, tmpdir):
        p = PfamDomtblout(domtblout_file)
        p.read()

        gff_with_extra = tmpdir.join("mixed.gff")
        gff_with_extra.write(
            "##gff-version 3\n"
            "chr1\tAUGUSTUS\tgene\t1000\t2000\t.\t+\t.\tID=gene1\n"
            "chr1\tAUGUSTUS\tmRNA\t1000\t2000\t.\t+\t.\tID=gene1.t1;Parent=gene1\n"
        )
        outfile = str(tmpdir.join("annotated.gff"))
        p.annotate_gff(str(gff_with_extra), outfile)

        with open(outfile) as f:
            content = f.read()
        assert "mRNA" in content
        assert "##gff-version 3" in content
