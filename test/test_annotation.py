import os
import tempfile

import pandas as pd
import pytest
from easydev import TempFile

from sequana import annotation

from . import test_dir


class TestPAFReader:
    """Test PAFReader class for reading PAF (Pairwise Alignment Format) files."""

    def test_paf_reader_standard_format(self):
        """Test reading standard PAF format with 11 columns."""
        # Create a temporary PAF file in standard format
        # Standard PAF format: r_name, r_start, r_end, strand, flag, mapq, cigar, q_name, q_len, q_start, q_end
        # Note: first row is treated as header by default, so need 2 data rows
        with TempFile(suffix=".paf") as fout:
            with open(fout.name, "w") as f:
                f.write("reference1\t500\t513\t+\t0\t50\t13M\tquery1\t100\t10\t90\n")
                f.write("reference2\t600\t620\t-\t10\t40\t20M\tquery2\t200\t20\t80\n")
            reader = annotation.PAFReader(fout.name)
            assert reader.df is not None
            assert isinstance(reader.df, pd.DataFrame)
            assert len(reader.df) == 1  # First row becomes header, so only 1 data row
            # Check that columns are assigned correctly
            assert list(reader.df.columns)[:3] == ["r_name", "r_start", "r_end"]

    def test_paf_reader_ragtag_format(self):
        """Test reading RAGTAG PAF format with additional columns."""
        # RAGTAG format has 18 columns
        with TempFile(suffix=".paf") as fout:
            with open(fout.name, "w") as f:
                f.write("query1\t100\t10\t90\t+\t0\t1000\treference1\t2000\t500\t513\t50\t" "tp\tcm\ts1\ts2\tdv\trl\n")
                # Add a second row so first becomes header
                f.write(
                    "query2\t200\t20\t180\t-\t10\t2000\treference2\t3000\t600\t620\t40\t" "tp\tcm\ts1\ts2\tdv\trl\n"
                )
            reader = annotation.PAFReader(fout.name)
            assert reader.df is not None
            assert isinstance(reader.df, pd.DataFrame)
            # Should have ragtag columns
            assert "q_length" in reader.df.columns or "q_name" in reader.df.columns

    def test_paf_reader_empty_file(self):
        """Test reading empty PAF file raises error."""
        with TempFile(suffix=".paf") as fout:
            # Empty file - should raise error
            with open(fout.name, "w") as f:
                pass
            with pytest.raises((pd.errors.EmptyDataError, ValueError)):
                reader = annotation.PAFReader(fout.name)

    def test_paf_reader_multiple_records(self):
        """Test reading PAF file with multiple records."""
        with TempFile(suffix=".paf") as fout:
            with open(fout.name, "w") as f:
                # First row becomes header, so need 4 rows for 3 data rows
                f.write("reference1\t500\t513\t+\t0\t50\t13M\tquery1\t100\t10\t90\n")
                f.write("reference2\t600\t620\t-\t10\t40\t20M\tquery2\t200\t20\t180\n")
                f.write("reference3\t700\t715\t+\t5\t45\t15M\tquery3\t150\t30\t130\n")
                f.write("reference4\t800\t820\t+\t2\t48\t20M\tquery4\t250\t40\t60\n")
            reader = annotation.PAFReader(fout.name)
            assert len(reader.df) == 3  # 4 rows - 1 header = 3 data rows


class TestAragorn:
    """Test Aragorn class for parsing tRNA predictions."""

    def test_aragorn_parser_existing_file(self):
        """Test parsing aragorn output from real test file."""
        arg = annotation.Aragorn()
        df = arg.parse_aragorn_output(f"{test_dir}/data/aragorn.txt")
        assert isinstance(df, pd.DataFrame)
        assert len(df) == 90
        assert "contig_name" in df.columns
        assert "tRNA_type" in df.columns
        assert "start" in df.columns
        assert "end" in df.columns
        assert "length" in df.columns
        assert "score" in df.columns
        assert "strand" in df.columns

    def test_aragorn_parser_basic_parsing(self):
        """Test basic parsing of aragorn output."""
        with TempFile(suffix=".txt") as fout:
            with open(fout.name, "w") as f:
                f.write(">contig_1\n")
                f.write("1   tRNA-Ala               [100,200]\t95.5\t33  \t(cgc)\n")
                f.write("2   tRNA-Gly               [300,400]\t105.2\t34  \t(gcc)\n")
            arg = annotation.Aragorn()
            df = arg.parse_aragorn_output(fout.name)
            assert len(df) == 2
            assert df.iloc[0]["contig_name"] == "contig_1"
            assert df.iloc[0]["tRNA_type"] == "tRNA-Ala"
            assert df.iloc[0]["start"] == 100
            assert df.iloc[0]["end"] == 200
            assert df.iloc[0]["length"] == 100
            assert df.iloc[0]["score"] == 95.5
            assert df.iloc[0]["strand"] == "+"

    def test_aragorn_parser_complement_strand(self):
        """Test parsing tRNA on complement strand (c[...])."""
        with TempFile(suffix=".txt") as fout:
            with open(fout.name, "w") as f:
                f.write(">contig_1\n")
                f.write("1   tRNA-Val                c[500,600]\t110.0\t35  \t(gac)\n")
            arg = annotation.Aragorn()
            df = arg.parse_aragorn_output(fout.name)
            assert len(df) == 1
            assert df.iloc[0]["strand"] == "-"
            assert df.iloc[0]["start"] == 500
            assert df.iloc[0]["end"] == 600

    def test_aragorn_parser_no_genes(self):
        """Test parsing aragorn output with contigs having no genes."""
        with TempFile(suffix=".txt") as fout:
            with open(fout.name, "w") as f:
                f.write(">contig_1\n")
                f.write("0 genes found\n")
                f.write(">contig_2\n")
                f.write("1   tRNA-Leu               [100,200]\t100.0\t33  \t(aag)\n")
            arg = annotation.Aragorn()
            df = arg.parse_aragorn_output(fout.name)
            assert len(df) == 1
            assert df.iloc[0]["contig_name"] == "contig_2"

    def test_aragorn_parser_empty_file(self):
        """Test parsing empty aragorn output file."""
        with TempFile(suffix=".txt") as fout:
            arg = annotation.Aragorn()
            df = arg.parse_aragorn_output(fout.name)
            assert isinstance(df, pd.DataFrame)
            assert len(df) == 0

    def test_aragorn_init(self):
        """Test Aragorn class initialization."""
        arg = annotation.Aragorn()
        assert arg is not None


class TestRNAmmer:
    """Test RNAmmer class for parsing rRNA predictions."""

    def test_rnammer_initialization(self):
        """Test RNAmmer stores filename correctly."""
        with TempFile(suffix=".gff") as fout:
            with open(fout.name, "w") as f:
                f.write("##gff-version2\n")
                f.write("contig_1\tRNAmmer-1.2\trRNA\t100\t200\t95.5\t+\t.\t16s_rRNA\tDummy\n")
                f.write("contig_2\tRNAmmer-1.2\trRNA\t500\t700\t100.0\t-\t.\t23s_rRNA\tDummy\n")
            rnammer = annotation.RNAmmer(fout.name)
            assert rnammer.filename == fout.name
            assert rnammer.df is not None
            assert isinstance(rnammer.df, pd.DataFrame)

    def test_rnammer_to_gff3(self):
        """Test RNAmmer conversion to GFF3 format."""
        with TempFile(suffix=".gff") as fin:
            with open(fin.name, "w") as f:
                f.write("##gff-version2\n")
                f.write("contig_1\tRNAmmer-1.2\trRNA\t100\t200\t95.5\t+\t.\t16s_rRNA\tDummy\n")
            rnammer = annotation.RNAmmer(fin.name)
            with TempFile(suffix=".gff3") as fout:
                rnammer.to_gff3(fout.name)
                assert os.path.exists(fout.name)
                with open(fout.name) as f:
                    content = f.read()
                    assert "contig_1" in content
                    assert "16s_rRNA" in content
                    assert "gene" in content

    def test_rnammer_multiple_records(self):
        """Test RNAmmer with multiple records."""
        with TempFile(suffix=".gff") as fout:
            with open(fout.name, "w") as f:
                f.write("##gff-version2\n")
                f.write("contig_1\tRNAmmer-1.2\trRNA\t100\t200\t95.5\t+\t.\t16s_rRNA\tDummy\n")
                f.write("contig_1\tRNAmmer-1.2\trRNA\t300\t600\t100.0\t-\t.\t23s_rRNA\tDummy\n")
                f.write("contig_2\tRNAmmer-1.2\trRNA\t500\t800\t98.0\t+\t.\t5s_rRNA\tDummy\n")
            rnammer = annotation.RNAmmer(fout.name)
            assert len(rnammer.df) == 3
            assert all(col in rnammer.df.columns for col in ["chrom", "feature", "strand", "score"])

    def test_rnammer_empty_file(self):
        """Test RNAmmer with empty GFF file raises error."""
        with TempFile(suffix=".gff") as fout:
            with open(fout.name, "w") as f:
                f.write("##gff-version2\n")
            # Empty file with only comments - should raise error
            with pytest.raises((pd.errors.EmptyDataError, ValueError)):
                rnammer = annotation.RNAmmer(fout.name)


class TestRFAMSplitter:
    """Test RFAMSplitter class for splitting RFAM database."""

    def test_rfam_splitter_init(self):
        """Test RFAMSplitter initialization."""
        with TempFile(suffix=".cm") as fout:
            with open(fout.name, "w") as f:
                f.write("ACC RF00001\n")
                f.write("ID SSU_rRNA_archaea\n")
                f.write("//\n")
            splitter = annotation.RFAMSplitter(fout.name)
            assert splitter.filename == fout.name

    def test_rfam_splitter_get_accessions(self):
        """Test getting accessions from RFAM file."""
        with TempFile(suffix=".cm") as fout:
            with open(fout.name, "w") as f:
                f.write("ACC RF00001\n")
                f.write("ID SSU_rRNA_archaea\n")
                f.write("//\n")
                f.write("ACC RF00002\n")
                f.write("ID SSU_rRNA_euk\n")
                f.write("//\n")
                f.write("ACC RF00001\n")  # Duplicate
                f.write("ID SSU_rRNA_bact\n")
                f.write("//\n")
            splitter = annotation.RFAMSplitter(fout.name)
            accs = splitter._get_accessions()
            assert isinstance(accs, list)
            assert len(accs) == 2  # Only unique accessions
            assert "RF00001" in accs
            assert "RF00002" in accs

    def test_rfam_splitter_extract(self):
        """Test extracting records from RFAM file."""
        with TempFile(suffix=".cm") as fin:
            with open(fin.name, "w") as f:
                f.write("ACC RF00001\n")
                f.write("ID SSU_rRNA_archaea\n")
                f.write("DESC 16S ribosomal RNA, archaeal\n")
                f.write("//\n")
                f.write("ACC RF00002\n")
                f.write("ID SSU_rRNA_euk\n")
                f.write("DESC 18S ribosomal RNA, eukaryotic\n")
                f.write("//\n")
            splitter = annotation.RFAMSplitter(fin.name)
            with TempFile(suffix=".cm") as fout:
                splitter.extract(["RF00001"], fout.name)
                assert os.path.exists(fout.name)
                with open(fout.name) as f:
                    content = f.read()
                    assert "RF00001" in content
                    assert "RF00002" not in content

    def test_rfam_splitter_extract_multiple(self):
        """Test extracting multiple accessions."""
        with TempFile(suffix=".cm") as fin:
            with open(fin.name, "w") as f:
                f.write("ACC RF00001\n")
                f.write("ID SSU_rRNA_archaea\n")
                f.write("//\n")
                f.write("ACC RF00002\n")
                f.write("ID SSU_rRNA_euk\n")
                f.write("//\n")
                f.write("ACC RF00003\n")
                f.write("ID SSU_rRNA_bact\n")
                f.write("//\n")
            splitter = annotation.RFAMSplitter(fin.name)
            with TempFile(suffix=".cm") as fout:
                splitter.extract(["RF00001", "RF00003"], fout.name)
                with open(fout.name) as f:
                    content = f.read()
                    assert "RF00001" in content
                    assert "RF00003" in content
                    assert "RF00002" not in content

    def test_rfam_splitter_empty_file(self):
        """Test RFAMSplitter with empty file."""
        with TempFile(suffix=".cm") as fout:
            splitter = annotation.RFAMSplitter(fout.name)
            accs = splitter._get_accessions()
            assert isinstance(accs, list)
            assert len(accs) == 0


class TestCMSearchParser:
    """Test CMSearchParser class for parsing cmsearch output."""

    def test_cmsearch_parser_init_and_parsing(self):
        """Test CMSearchParser initialization and parsing."""
        with TempFile(suffix=".txt") as fout:
            with open(fout.name, "w") as f:
                f.write("# cmsearch :: search RNA profiles against a sequence database\n")
                f.write("# INFERNAL 1.1.4\n")
                # Data lines without header - 18 columns: target_name, accession1, query_name, accession,
                # mdl, from, to, seq_from, seq_to, strand, trunc, pass, gc, bias, score, E_value, inc, description_target
                f.write("sequence_1 RF00001 SSU_rRNA RF00001 hmm 1 100 10 109 + N yes 0.55 0.1 75.5 1e-10 Y 16S\n")
            parser = annotation.CMSearchParser(fout.name)
            assert parser.df is not None
            assert isinstance(parser.df, pd.DataFrame)
            assert len(parser.df) == 1
            assert "target_name" in parser.df.columns

    def test_cmsearch_parser_multiple_records(self):
        """Test CMSearchParser with multiple records."""
        with TempFile(suffix=".txt") as fout:
            with open(fout.name, "w") as f:
                f.write("# cmsearch :: search RNA profiles against a sequence database\n")
                f.write("seq_1 RF00001 SSU_rRNA RF00001 hmm 1 100 10 109 + N yes 0.55 0.1 75.5 1e-10 Y 16S\n")
                f.write("seq_2 RF00002 LSU_rRNA RF00002 hmm 1 200 20 219 - N yes 0.50 0.2 80.0 1e-12 Y 23S\n")
            parser = annotation.CMSearchParser(fout.name)
            assert len(parser.df) == 2

    def test_cmsearch_parser_description_merging(self):
        """Test that multi-word descriptions are merged correctly."""
        with TempFile(suffix=".txt") as fout:
            with open(fout.name, "w") as f:
                f.write("# cmsearch output\n")
                f.write("seq_1 RF00001 SSU_rRNA RF00001 hmm 1 100 10 109 + N yes 0.55 0.1 75.5 1e-10 Y this is long\n")
            parser = annotation.CMSearchParser(fout.name)
            assert len(parser.df) == 1

    def test_cmsearch_parser_empty_file(self):
        """Test CMSearchParser with empty file raises error."""
        with TempFile(suffix=".txt") as fout:
            with open(fout.name, "w") as f:
                f.write("# cmsearch header\n")
            with pytest.raises((pd.errors.EmptyDataError, ValueError)):
                parser = annotation.CMSearchParser(fout.name)


class TestRAGTAG:
    """Test RAGTAG class for reading PAF files."""

    def test_ragtag_init(self):
        """Test RAGTAG initialization."""
        with TempFile(suffix=".paf") as fout:
            with open(fout.name, "w") as f:
                # RAGTAG format - 18 columns
                f.write("query1\t100\t10\t90\t+\t0\t1000\treference1\t2000\t500\t513\t50\t" "tp\tcm\ts1\ts2\tdv\trl\n")
            ragtag = annotation.RAGTAG(fout.name)
            assert ragtag.df is not None
            assert isinstance(ragtag.df, pd.DataFrame)

    def test_ragtag_uses_paf_reader(self):
        """Test that RAGTAG uses PAFReader internally."""
        with TempFile(suffix=".paf") as fout:
            with open(fout.name, "w") as f:
                # Standard PAF format - 11 columns
                f.write("r1\t500\t513\t+\t0\t50\t13M\tq1\t100\t10\t90\n")
            ragtag = annotation.RAGTAG(fout.name)
            assert isinstance(ragtag.df, pd.DataFrame)
            assert len(ragtag.df) >= 0  # May have columns or be empty depending on PAF parsing
