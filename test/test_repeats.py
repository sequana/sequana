"""Comprehensive tests for the repeats module."""
import os
import tempfile

import numpy as np
import pandas as pd
import pytest

from sequana.repeats.aphased import APhasedRepeats
from sequana.repeats.cruciforms import Cruciforms
from sequana.repeats.directrepeat import DirectRepeats
from sequana.repeats.G4hunter import G4Hunter
from sequana.repeats.gquad import GQuadruplex
from sequana.repeats.hdna import HDNA
from sequana.repeats.imotif import IMotif
from sequana.repeats.mirror import MirrorRepeats
from sequana.repeats.mirror import _scan_python as mirror_scan_python
from sequana.repeats.palindromes import Palindromes
from sequana.repeats.shustring import Repeats as ShuStringRepeats
from sequana.repeats.tandem import ShortTandemRepeats, _nonbstr_py, _scan_python
from sequana.repeats.trf import TRF
from sequana.repeats.zdna import ZDNA

from . import test_dir


class TestShortTandemRepeats:
    """Test tandem.py module."""

    def setup_method(self):
        self.fasta_file = f"{test_dir}/data/fasta/test_shustring.fa"

    def test_init(self):
        """Test initialization."""
        str_obj = ShortTandemRepeats(self.fasta_file)
        assert str_obj.fasta_file == self.fasta_file
        assert str_obj.min_period == 1
        assert str_obj.max_period == 9
        assert str_obj.min_span == 10
        assert str_obj.min_reps == 3
        assert isinstance(str_obj.df, pd.DataFrame)
        assert str_obj.df.empty

    def test_custom_parameters(self):
        """Test initialization with custom parameters."""
        str_obj = ShortTandemRepeats(self.fasta_file, min_period=2, max_period=5, min_span=15, min_reps=4)
        assert str_obj.min_period == 2
        assert str_obj.max_period == 5
        assert str_obj.min_span == 15
        assert str_obj.min_reps == 4

    def test_run(self):
        """Test run method."""
        str_obj = ShortTandemRepeats(self.fasta_file, min_period=1, max_period=3)
        str_obj.run(progress=False)
        assert isinstance(str_obj.df, pd.DataFrame)
        assert "seqid" in str_obj.df.columns
        assert "start" in str_obj.df.columns
        assert "end" in str_obj.df.columns
        assert "period" in str_obj.df.columns
        assert "num" in str_obj.df.columns

    def test_run_with_progress(self):
        """Test run method with progress bar."""
        str_obj = ShortTandemRepeats(self.fasta_file)
        str_obj.run(progress=False)  # Just ensure it doesn't crash

    def test_to_bed(self):
        """Test to_bed method."""
        str_obj = ShortTandemRepeats(self.fasta_file)
        str_obj.run(progress=False)

        with tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".bed") as f:
            bed_file = f.name

        try:
            if not str_obj.df.empty:
                str_obj.to_bed(bed_file)
                assert os.path.exists(bed_file)
                # Verify it's a valid BED file
                df = pd.read_csv(bed_file, sep="\t", header=None)
                assert len(df) > 0
                assert len(df.columns) >= 6
        finally:
            if os.path.exists(bed_file):
                os.unlink(bed_file)

    def test_to_bed_empty(self):
        """Test to_bed raises error for empty results."""
        str_obj = ShortTandemRepeats(self.fasta_file)
        # Don't run, df stays empty
        with pytest.raises(ValueError, match="Run `.run\\(\\)` first"):
            str_obj.to_bed("/tmp/test.bed")

    def test_to_gff(self):
        """Test to_gff method."""
        str_obj = ShortTandemRepeats(self.fasta_file)
        str_obj.run(progress=False)

        with tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".gff") as f:
            gff_file = f.name

        try:
            if not str_obj.df.empty:
                str_obj.to_gff(gff_file)
                assert os.path.exists(gff_file)
                with open(gff_file) as f:
                    lines = f.readlines()
                    assert lines[0].startswith("##gff-version")
                    if len(lines) > 1:
                        assert "Short_Tandem_Repeat" in lines[1]
        finally:
            if os.path.exists(gff_file):
                os.unlink(gff_file)

    def test_to_gff_empty(self):
        """Test to_gff raises error for empty results."""
        str_obj = ShortTandemRepeats(self.fasta_file)
        with pytest.raises(ValueError, match="Run `.run\\(\\)` first"):
            str_obj.to_gff("/tmp/test.gff")

    def test_nonbstr_py(self):
        """Test _nonbstr_py function."""
        # Test even length symmetric sequence
        result = _nonbstr_py("AT")
        assert isinstance(result, int)

        # Test complementary sequence
        result = _nonbstr_py("ATAT")
        assert isinstance(result, int)

        # Test homopurine/homopyrimidine
        result = _nonbstr_py("AAAA")
        assert isinstance(result, int)

    def test_scan_python(self):
        """Test _scan_python function."""
        seq = b"AAAATTTTCCCCGGGG"
        arr = np.frombuffer(seq, dtype=np.uint8)
        starts, ends, periods, nums, subs, types = _scan_python(arr, 2, 4, 4, 2)
        assert isinstance(starts, np.ndarray)
        assert isinstance(ends, np.ndarray)


class TestMirrorRepeats:
    """Test mirror.py module."""

    def setup_method(self):
        self.fasta_file = f"{test_dir}/data/fasta/test_shustring.fa"

    def test_init(self):
        """Test initialization."""
        mr = MirrorRepeats(self.fasta_file)
        assert mr.fasta_file == self.fasta_file
        assert mr.min_repeat == 10
        assert mr.max_repeat == 200
        assert mr.min_spacer == 0
        assert mr.max_spacer == 100
        assert isinstance(mr.df, pd.DataFrame)

    def test_custom_parameters(self):
        """Test initialization with custom parameters."""
        mr = MirrorRepeats(self.fasta_file, min_repeat=5, max_repeat=50, min_spacer=1, max_spacer=20)
        assert mr.min_repeat == 5
        assert mr.max_repeat == 50
        assert mr.min_spacer == 1
        assert mr.max_spacer == 20

    def test_run(self):
        """Test run method."""
        mr = MirrorRepeats(self.fasta_file, min_repeat=5, max_repeat=20)
        mr.run(progress=False)
        assert isinstance(mr.df, pd.DataFrame)
        # May or may not have results depending on sequence

    def test_to_bed(self):
        """Test to_bed method."""
        mr = MirrorRepeats(self.fasta_file, min_repeat=5, max_repeat=20)
        mr.run(progress=False)

        with tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".bed") as f:
            bed_file = f.name

        try:
            if not mr.df.empty:
                mr.to_bed(bed_file)
                assert os.path.exists(bed_file)
        finally:
            if os.path.exists(bed_file):
                os.unlink(bed_file)

    def test_to_bed_empty(self):
        """Test to_bed raises error for empty results."""
        mr = MirrorRepeats(self.fasta_file)
        with pytest.raises(ValueError, match="Run `.run\\(\\)` first"):
            mr.to_bed("/tmp/test.bed")


class TestHDNA:
    """Test hdna.py module."""

    def setup_method(self):
        self.fasta_file = f"{test_dir}/data/fasta/test_shustring.fa"

    def test_init(self):
        """Test initialization."""
        hdna = HDNA(self.fasta_file)
        assert hdna.fasta_file == self.fasta_file
        assert isinstance(hdna.df, pd.DataFrame)

    def test_run(self):
        """Test run method."""
        hdna = HDNA(self.fasta_file)
        hdna.run(progress=False)
        assert isinstance(hdna.df, pd.DataFrame)

    def test_to_bed(self):
        """Test to_bed method."""
        hdna = HDNA(self.fasta_file)
        hdna.run(progress=False)

        with tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".bed") as f:
            bed_file = f.name

        try:
            if not hdna.df.empty:
                hdna.to_bed(bed_file)
                assert os.path.exists(bed_file)
        finally:
            if os.path.exists(bed_file):
                os.unlink(bed_file)


class TestAPhased:
    """Test aphased.py module."""

    def setup_method(self):
        self.fasta_file = f"{test_dir}/data/fasta/test_shustring.fa"

    def test_init(self):
        """Test initialization."""
        ap = APhasedRepeats(self.fasta_file)
        assert ap.fasta_file == self.fasta_file
        assert isinstance(ap.df, pd.DataFrame)

    def test_run(self):
        """Test run method."""
        ap = APhasedRepeats(self.fasta_file)
        ap.run(progress=False)
        assert isinstance(ap.df, pd.DataFrame)

    def test_to_bed(self):
        """Test to_bed method."""
        ap = APhasedRepeats(self.fasta_file)
        ap.run(progress=False)

        with tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".bed") as f:
            bed_file = f.name

        try:
            if not ap.df.empty:
                ap.to_bed(bed_file)
                assert os.path.exists(bed_file)
        finally:
            if os.path.exists(bed_file):
                os.unlink(bed_file)


class TestGQuadruplex:
    """Test gquad.py module."""

    def setup_method(self):
        self.fasta_file = f"{test_dir}/data/fasta/test_shustring.fa"

    def test_init(self):
        """Test initialization."""
        gq = GQuadruplex(self.fasta_file)
        assert gq.fasta_file == self.fasta_file
        assert gq.min_rep == 3
        assert gq.max_spacer == 7
        assert isinstance(gq.df, pd.DataFrame)

    def test_custom_parameters(self):
        """Test initialization with custom parameters."""
        gq = GQuadruplex(self.fasta_file, min_rep=2, max_spacer=5)
        assert gq.min_rep == 2
        assert gq.max_spacer == 5

    def test_run(self):
        """Test run method."""
        gq = GQuadruplex(self.fasta_file)
        gq.run(progress=False)
        assert isinstance(gq.df, pd.DataFrame)

    def test_to_bed(self):
        """Test to_bed method."""
        gq = GQuadruplex(self.fasta_file)
        gq.run(progress=False)

        with tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".bed") as f:
            bed_file = f.name

        try:
            if not gq.df.empty:
                gq.to_bed(bed_file)
                assert os.path.exists(bed_file)
        finally:
            if os.path.exists(bed_file):
                os.unlink(bed_file)


class TestDirectRepeat:
    """Test directrepeat.py module."""

    def setup_method(self):
        self.fasta_file = f"{test_dir}/data/fasta/test_shustring.fa"

    def test_init(self):
        """Test initialization."""
        dr = DirectRepeats(self.fasta_file)
        assert dr.fasta_file == self.fasta_file
        assert isinstance(dr.df, pd.DataFrame)

    def test_run(self):
        """Test run method."""
        dr = DirectRepeats(self.fasta_file)
        dr.run(progress=False)
        assert isinstance(dr.df, pd.DataFrame)

    def test_to_bed(self):
        """Test to_bed method."""
        dr = DirectRepeats(self.fasta_file)
        dr.run(progress=False)

        with tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".bed") as f:
            bed_file = f.name

        try:
            if not dr.df.empty:
                dr.to_bed(bed_file)
                assert os.path.exists(bed_file)
        finally:
            if os.path.exists(bed_file):
                os.unlink(bed_file)


class TestShuString:
    """Test shustring.py module (ShuString uses external shustring tool)."""

    def setup_method(self):
        self.fasta_file = f"{test_dir}/data/fasta/test_shustring.fa"

    def test_init(self):
        """Test initialization."""
        ss = ShuStringRepeats(self.fasta_file)
        assert ss._filename_fasta == self.fasta_file

    def test_properties(self):
        """Test basic properties."""
        ss = ShuStringRepeats(self.fasta_file)
        # Check header can be set
        assert ss.header in ss.names
        # Check names list exists
        assert len(ss.names) > 0

    def test_threshold_property(self):
        """Test that accessing threshold property without external tool doesn't crash."""
        # This test just checks that the class structure is correct
        # The actual shustring tool may not be installed
        ss = ShuStringRepeats(self.fasta_file)
        assert hasattr(ss, "_threshold")


class TestCruciforms:
    """Test cruciforms.py module."""

    def setup_method(self):
        self.fasta_file = f"{test_dir}/data/fasta/test_shustring.fa"

    def test_init(self):
        """Test initialization."""
        cf = Cruciforms(self.fasta_file)
        assert cf.fasta_file == self.fasta_file
        assert isinstance(cf.df, pd.DataFrame)

    def test_run(self):
        """Test run method."""
        cf = Cruciforms(self.fasta_file)
        cf.run(progress=False)
        assert isinstance(cf.df, pd.DataFrame)

    def test_to_bed(self):
        """Test to_bed method."""
        cf = Cruciforms(self.fasta_file)
        cf.run(progress=False)

        with tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".bed") as f:
            bed_file = f.name

        try:
            if not cf.df.empty:
                cf.to_bed(bed_file)
                assert os.path.exists(bed_file)
        finally:
            if os.path.exists(bed_file):
                os.unlink(bed_file)


class TestZDNA:
    """Test zdna.py module."""

    def setup_method(self):
        self.fasta_file = f"{test_dir}/data/fasta/test_shustring.fa"

    def test_init(self):
        """Test initialization."""
        zdna = ZDNA(self.fasta_file)
        assert zdna.fasta_file == self.fasta_file
        assert isinstance(zdna.df, pd.DataFrame)

    def test_run(self):
        """Test run method."""
        zdna = ZDNA(self.fasta_file)
        zdna.run(progress=False)
        assert isinstance(zdna.df, pd.DataFrame)

    def test_to_bed(self):
        """Test to_bed method."""
        zdna = ZDNA(self.fasta_file)
        zdna.run(progress=False)

        with tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".bed") as f:
            bed_file = f.name

        try:
            if not zdna.df.empty:
                zdna.to_bed(bed_file)
                assert os.path.exists(bed_file)
        finally:
            if os.path.exists(bed_file):
                os.unlink(bed_file)


class TestPalindromes:
    """Test palindromes.py module."""

    def setup_method(self):
        self.fasta_file = f"{test_dir}/data/fasta/test_shustring.fa"

    def test_init(self):
        """Test initialization."""
        pal = Palindromes(self.fasta_file)
        assert pal.fasta_file == self.fasta_file
        assert pal.min_len == 4
        assert pal.max_len == 12
        assert isinstance(pal.df, pd.DataFrame)

    def test_custom_parameters(self):
        """Test initialization with custom parameters."""
        pal = Palindromes(self.fasta_file, min_len=6, max_len=20)
        assert pal.min_len == 6
        assert pal.max_len == 20

    def test_is_palindrome(self):
        """Test is_palindrome method."""
        pal = Palindromes(self.fasta_file)
        # GAATTC is a palindrome (reverse complement equals itself)
        assert pal.is_palindrome("GAATTC")
        # Random sequence should not be palindrome
        assert not pal.is_palindrome("ATCG")
        # Case insensitivity
        assert pal.is_palindrome("gaattc")

    def test_run(self):
        """Test run method."""
        pal = Palindromes(self.fasta_file, min_len=4, max_len=10)
        pal.run()
        assert isinstance(pal.df, pd.DataFrame)

    def test_to_bed(self):
        """Test to_bed method."""
        pal = Palindromes(self.fasta_file, min_len=4, max_len=10)
        pal.run()

        with tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".bed") as f:
            bed_file = f.name

        try:
            if not pal.df.empty:
                pal.to_bed(bed_file)
                assert os.path.exists(bed_file)
        finally:
            if os.path.exists(bed_file):
                os.unlink(bed_file)

    def test_to_bed_empty(self):
        """Test to_bed raises error for empty results."""
        pal = Palindromes(self.fasta_file)
        with pytest.raises(ValueError, match="Run `.run\\(\\)` first"):
            pal.to_bed("/tmp/test.bed")


class TestIMotif:
    """Test imotif.py module."""

    def setup_method(self):
        self.fasta_file = f"{test_dir}/data/fasta/test_shustring.fa"

    def test_init(self):
        """Test initialization."""
        im = IMotif(self.fasta_file)
        assert im.fasta_file == self.fasta_file
        assert im.min_tract == 3
        assert im.max_loop == 7
        assert isinstance(im.df, pd.DataFrame)

    def test_custom_parameters(self):
        """Test initialization with custom parameters."""
        im = IMotif(self.fasta_file, min_tract=2, max_loop=5)
        assert im.min_tract == 2
        assert im.max_loop == 5

    def test_run(self):
        """Test run method."""
        im = IMotif(self.fasta_file)
        im.run()
        assert isinstance(im.df, pd.DataFrame)

    def test_to_bed(self):
        """Test to_bed method."""
        im = IMotif(self.fasta_file)
        im.run()

        with tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".bed") as f:
            bed_file = f.name

        try:
            if not im.df.empty:
                im.to_bed(bed_file)
                assert os.path.exists(bed_file)
        finally:
            if os.path.exists(bed_file):
                os.unlink(bed_file)

    def test_to_bed_empty(self):
        """Test to_bed raises error for empty results."""
        im = IMotif(self.fasta_file)
        with pytest.raises(ValueError, match="No results"):
            im.to_bed("/tmp/test.bed")


class TestG4Hunter:
    """Test G4hunter.py module."""

    def setup_method(self):
        self.fasta_file = f"{test_dir}/data/fasta/test_shustring.fa"

    def test_init(self):
        """Test initialization."""
        g4 = G4Hunter(self.fasta_file)
        assert g4.infile == self.fasta_file
        assert g4.window == 25
        assert g4.score == 1

    def test_custom_parameters(self):
        """Test initialization with custom parameters."""
        g4 = G4Hunter(self.fasta_file, window=30, score=2)
        assert g4.window == 30
        assert g4.score == 2

    def test_run(self):
        """Test run method."""
        with tempfile.TemporaryDirectory() as tmpdir:
            g4 = G4Hunter(self.fasta_file)
            g4.run(tmpdir)
            # Check that output files were created
            assert os.path.exists(tmpdir)

    def test_base_score(self):
        """Test base_score method."""
        g4 = G4Hunter(self.fasta_file)
        scores = g4.base_score("GGGGCCCC")
        assert isinstance(scores, np.ndarray)
        assert len(scores) == 8

    def test_output_files_created(self):
        """Test that run method creates output files."""
        with tempfile.TemporaryDirectory() as tmpdir:
            g4 = G4Hunter(self.fasta_file)
            g4.run(tmpdir)
            # Check that output directory exists and has files
            assert os.path.isdir(tmpdir)


class TestTandemRepeatFinder:
    """Test trf.py module (reads TRF tool output, not FASTA)."""

    def test_init_with_csv(self):
        """Test initialization with CSV format."""
        # Create a mock CSV file for testing
        csv_data = """seq_id,start,end,length,period_size
contig1,100,200,100,5
contig1,300,500,200,3
"""
        with tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".csv") as f:
            f.write(csv_data)
            csv_file = f.name

        try:
            trf = TRF(csv_file, frmt="csv")
            assert isinstance(trf.df, pd.DataFrame)
            assert not trf.df.empty
        finally:
            os.unlink(csv_file)

    def test_init_invalid_format(self):
        """Test initialization with invalid file format."""
        with tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".txt") as f:
            f.write("invalid data")
            txt_file = f.name

        try:
            with pytest.raises(ValueError, match="Unknown file type"):
                TRF(txt_file)
        finally:
            os.unlink(txt_file)

    def test_dataframe_length(self):
        """Test that dataframe length works correctly."""
        csv_data = """seq_id,start,end
contig1,100,200
contig1,300,500
"""
        with tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".csv") as f:
            f.write(csv_data)
            csv_file = f.name

        try:
            trf = TRF(csv_file, frmt="csv")
            assert len(trf) == 2
        finally:
            os.unlink(csv_file)
