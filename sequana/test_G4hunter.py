import os
import tempfile
from pathlib import Path

import pytest

from sequana.G4hunter import G4HunterReader


class TestG4HunterReaderLoadMergedData:
    def test_load_merged_data_basic(self):
        """Test loading a basic merged data file with one sequence"""
        with tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".txt") as f:
            f.write(">seq1\n")
            f.write("Start\tEnd\tSequence\tLength\tScore\n")
            f.write("0\t20\tACGTACGTACGTACGTACGT\t20\t1.5\n")
            f.write("25\t45\tTGCATGCATGCATGCATGCA\t20\t-2.3\n")
            f.flush()

            reader = G4HunterReader()
            reader.load_merged_data(f.name)

            assert "seq1" in reader.data_merged
            assert len(reader.data_merged["seq1"]) == 2
            assert reader.data_merged["seq1"].iloc[0]["Start"] == 0
            assert reader.data_merged["seq1"].iloc[1]["Score"] == -2.3

        os.unlink(f.name)

    def test_load_merged_data_with_nbr_column(self):
        """Test loading merged data with NBR column (6 items per line)"""
        with tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".txt") as f:
            f.write(">seq1\n")
            f.write("Start\tEnd\tSequence\tLength\tScore\tNBR\n")
            f.write("0\t20\tACGT\t4\t2.5\t3\n")
            f.flush()

            reader = G4HunterReader()
            reader.load_merged_data(f.name)

            assert reader.data_merged["seq1"].iloc[0]["NBR"] == 3
            assert reader.data_merged["seq1"].iloc[0]["Score"] == 2.5

        os.unlink(f.name)

    def test_load_merged_data_multiple_sequences(self):
        """Test loading file with multiple sequences"""
        with tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".txt") as f:
            f.write(">seq1\n")
            f.write("Start\tEnd\tSequence\tLength\tScore\n")
            f.write("0\t20\tACGT\t4\t1.0\n")
            f.write(">seq2\n")
            f.write("Start\tEnd\tSequence\tLength\tScore\n")
            f.write("10\t30\tTGCA\t4\t-1.5\n")
            f.flush()

            reader = G4HunterReader()
            reader.load_merged_data(f.name)

            assert len(reader.data_merged) == 2
            assert "seq1" in reader.data_merged
            assert "seq2" in reader.data_merged
            assert len(reader.data_merged["seq1"]) == 1
            assert len(reader.data_merged["seq2"]) == 1

        os.unlink(f.name)

    def test_load_merged_data_empty_file(self):
        """Test loading an empty file"""
        with tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".txt") as f:
            f.flush()

            reader = G4HunterReader()
            reader.load_merged_data(f.name)

            assert len(reader.data_merged) == 0

        os.unlink(f.name)

    def test_load_merged_data_invalid_items_count(self):
        """Test that lines with invalid item count are skipped"""
        with tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".txt") as f:
            f.write(">seq1\n")
            f.write("Start\tEnd\tSequence\tLength\tScore\n")
            f.write("0\t20\tACGT\t4\n")  # Only 4 items, should be skipped
            f.write("10\t30\tTGCA\t4\t1.5\n")  # Valid line
            f.flush()

            reader = G4HunterReader()
            reader.load_merged_data(f.name)

            assert len(reader.data_merged["seq1"]) == 1

        os.unlink(f.name)

    def test_load_merged_data_whitespace_handling(self):
        """Test that whitespace around values is properly stripped"""
        with tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".txt") as f:
            f.write(">seq1\n")
            f.write("Start\tEnd\tSequence\tLength\tScore\n")
            f.write("  0  \t  20  \t  ACGT  \t  4  \t  1.5  \n")
            f.flush()

            reader = G4HunterReader()
            reader.load_merged_data(f.name)

            assert reader.data_merged["seq1"].iloc[0]["Start"] == 0
            assert reader.data_merged["seq1"].iloc[0]["Sequence"] == "ACGT"

        os.unlink(f.name)

    def test_load_merged_data_negative_scores(self):
        """Test handling of negative score values"""
        with tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".txt") as f:
            f.write(">seq1\n")
            f.write("Start\tEnd\tSequence\tLength\tScore\n")
            f.write("0\t20\tACGT\t4\t-3.2\n")
            f.flush()

            reader = G4HunterReader()
            reader.load_merged_data(f.name)

            assert reader.data_merged["seq1"].iloc[0]["Score"] == -3.2

        os.unlink(f.name)

    def test_load_merged_data_header_extraction(self):
        """Test that header ID is correctly extracted"""
        with tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".txt") as f:
            f.write(">chromosome_1 extra info\n")
            f.write("Start\tEnd\tSequence\tLength\tScore\n")
            f.write("0\t20\tACGT\t4\t1.0\n")
            f.flush()

            reader = G4HunterReader()
            reader.load_merged_data(f.name)

            assert "chromosome_1" in reader.data_merged

        os.unlink(f.name)
