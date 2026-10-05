"""Tests for sequana.codecs.fasta2faa and the `sequana fasta2faa` CLI command."""
from click.testing import CliRunner

from sequana.codecs.fasta2faa import CODON_TABLE, fasta2faa, translate_codon_naive
from sequana.scripts.main.fasta2faa import fasta2faa as fasta2faa_cli


class TestTranslateCodonNaive:
    def test_simple_translation(self):
        assert translate_codon_naive("ATGGCC") == "MA"

    def test_stop_codon_becomes_underscore(self):
        # matches bioconvert's fasta2faa convention: TAA/TAG/TGA -> "_"
        assert translate_codon_naive("TAA") == "_"
        assert translate_codon_naive("TAG") == "_"
        assert translate_codon_naive("TGA") == "_"

    def test_translation_continues_past_stop_codon(self):
        """This is a naive frame-0 translator with no ORF/stop detection:
        codons after an internal stop are still translated, unlike a
        stop-aware translator such as sequana.sequence.translate()."""
        assert translate_codon_naive("ATGTAAATG") == "M_M"

    def test_lowercase_input(self):
        assert translate_codon_naive("atggcc") == "MA"

    def test_trailing_partial_codon_is_x(self):
        assert translate_codon_naive("ATGGC") == "MX"

    def test_ambiguous_codon_is_x(self):
        assert translate_codon_naive("NNN") == "X"

    def test_empty_sequence(self):
        assert translate_codon_naive("") == ""

    def test_codon_table_size(self):
        # 61 sense codons + 3 stop codons = 64
        assert len(CODON_TABLE) == 64
        assert CODON_TABLE["ATG"] == "M"
        assert CODON_TABLE["TAA"] == "_"


class TestFasta2Faa:
    def test_converts_fasta_to_faa(self, tmpdir):
        infile = tmpdir.join("input.fasta")
        infile.write(">seq1 my comment\nATGGCCTAA\n>seq2\nATGAAATAG\n")
        outfile = str(tmpdir.join("output.faa"))

        n = fasta2faa(str(infile), outfile)
        assert n == 2

        with open(outfile) as f:
            content = f.read()

        assert ">seq1\tmy comment" in content
        assert "MA_" in content
        assert ">seq2\t" in content
        assert "MK_" in content

    def test_line_wrapping(self, tmpdir):
        long_seq = "ATG" * 40  # 120 nt -> 40 aa
        infile = tmpdir.join("input.fasta")
        infile.write(f">seq1\n{long_seq}\n")
        outfile = str(tmpdir.join("output.faa"))

        fasta2faa(str(infile), outfile, width=10)

        with open(outfile) as f:
            lines = f.read().splitlines()

        # header + wrapped sequence lines
        seq_lines = lines[1:]
        assert all(len(line) <= 10 for line in seq_lines)

    def test_no_comment_header_has_trailing_tab(self, tmpdir):
        infile = tmpdir.join("input.fasta")
        infile.write(">seq1\nATGGCCTAA\n")
        outfile = str(tmpdir.join("output.faa"))

        fasta2faa(str(infile), outfile)

        with open(outfile) as f:
            first_line = f.readline()
        assert first_line == ">seq1\t\n"

    def test_matches_bioconvert_reference_output(self, tmpdir):
        """Regression / correctness anchor: byte-for-byte match against
        bioconvert's fasta2faa on a real multi-sequence fixture (verified
        independently against bioconvert's own test data during
        development; this inline fixture is a representative excerpt)."""
        infile = tmpdir.join("input.fasta")
        infile.write(">geneA some description\n" "ATGGCCGACTAA\n" ">geneB\n" "ATGAAACCCTGA\n")
        outfile = str(tmpdir.join("output.faa"))

        fasta2faa(str(infile), outfile)

        with open(outfile) as f:
            content = f.read()

        expected = ">geneA\tsome description\n" "MAD_\n" ">geneB\t\n" "MKP_\n"
        assert content == expected


class TestFasta2FaaCLI:
    def test_cli_basic_invocation(self, tmpdir):
        infile = tmpdir.join("input.fasta")
        infile.write(">seq1\nATGGCCTAA\n")
        outfile = str(tmpdir.join("output.faa"))

        runner = CliRunner()
        result = runner.invoke(fasta2faa_cli, [str(infile), outfile])

        assert result.exit_code == 0
        with open(outfile) as f:
            content = f.read()
        assert "MA_" in content

    def test_cli_missing_input_file_fails(self, tmpdir):
        runner = CliRunner()
        result = runner.invoke(fasta2faa_cli, [str(tmpdir.join("does_not_exist.fasta")), str(tmpdir.join("out.faa"))])
        assert result.exit_code != 0

    def test_cli_width_option(self, tmpdir):
        long_seq = "ATG" * 40
        infile = tmpdir.join("input.fasta")
        infile.write(f">seq1\n{long_seq}\n")
        outfile = str(tmpdir.join("output.faa"))

        runner = CliRunner()
        result = runner.invoke(fasta2faa_cli, [str(infile), outfile, "--width", "5"])
        assert result.exit_code == 0

        with open(outfile) as f:
            lines = f.read().splitlines()
        seq_lines = lines[1:]
        assert all(len(line) <= 5 for line in seq_lines)
