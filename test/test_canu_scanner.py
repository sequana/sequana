"""Tests for sequana.canu_scanner."""
import matplotlib
import pytest

matplotlib.use("Agg")

from sequana.canu_scanner import CanuScanner


@pytest.fixture
def canu_dir(tmpdir):
    """Build a minimal synthetic Canu output directory tree."""
    base = tmpdir
    corr_gkp = base.mkdir("correction").mkdir("sample.gkpStore")
    mercounts = base.join("correction").mkdir("0-mercounts")
    corr2 = base.join("correction").mkdir("2-correction")
    trim_gkp = base.mkdir("trimming").mkdir("sample.gkpStore")

    load_dat = "header line 0\nheader line 1\nReads 100 500000 5 2000\n"
    corr_gkp.join("load.dat").write(load_dat)
    trim_gkp.join("load.dat").write(load_dat)

    reads_txt = "read1\tx\t1000\ty\tz\nread2\tx\t2000\ty\tz\nread3\tx\t1500\ty\tz\n"
    corr_gkp.join("reads.txt").write(reads_txt)
    trim_gkp.join("reads.txt").write(reads_txt)

    mercounts.join("sample.ms16.histogram").write(
        "1\t500\t0.1\t0.05\n2\t300\t0.3\t0.2\n3\t100\t0.6\t0.5\n4\t50\t1.0\t1.0\n"
    )

    corr2.join("sample.globalScores.stats").write("some overlap filtering stats text\n")
    corr2.join("sample.correction.summary").write("some correction summary text\n")
    corr2.join("sample.original-expected-corrected-length.dat").write(
        "read1\t1000\t1010\t1005\nread2\t2000\t1990\t1995\n"
    )

    return str(base)


class TestCanuScanner:
    def test_init_defaults(self):
        c = CanuScanner()
        assert c.data["tool"] == "sequana"
        assert c.data["module"] == "canu_scanner"
        assert c.data["correction"] == {}

    def test_getfile_raises_when_no_match(self, canu_dir):
        c = CanuScanner(canu_dir)
        with pytest.raises(FileNotFoundError):
            c.getfile("does_not_exist/*.foo")

    def test_getfile_raises_when_multiple_matches(self, canu_dir):
        c = CanuScanner(canu_dir)
        # both correction and trimming gkpStores have a load.dat
        with pytest.raises(FileNotFoundError):
            c.getfile("*/*.gkpStore/load.dat")

    def test_scan_correction(self, canu_dir):
        c = CanuScanner(canu_dir)
        c.scan_correction()
        assert c.data["correction"]["readsLoaded"] == {"reads": 100, "bp": 500000}
        assert c.data["correction"]["readsSkipped"] == {"reads": 5, "bp": 2000}

    def test_scan_trimming(self, canu_dir):
        c = CanuScanner(canu_dir)
        c.scan_trimming()
        assert c.data["trimming"]["readsLoaded"] == {"reads": 100, "bp": 500000}
        assert c.data["trimming"]["readsSkipped"] == {"reads": 5, "bp": 2000}

    def test_hist_read_length(self, canu_dir):
        c = CanuScanner(canu_dir)
        df = c.hist_read_length()
        assert len(df) == 3
        assert list(df.columns) == ["ID", 1, "read_length", 3, 4]

    def test_hist_trimming_read_length(self, canu_dir):
        c = CanuScanner(canu_dir)
        df = c.hist_trimming_read_length()
        assert len(df) == 3

    def test_plot_kmer_populates_data(self, canu_dir):
        c = CanuScanner(canu_dir)
        df = c.plot_kmer()
        assert len(df) == 4
        assert c.data["correction"]["largest mercount"] == 4
        assert c.data["correction"]["unique mers"] == 500
        assert c.data["correction"]["distinc mers"] == 950
        # the previous implementation stored this under an empty-string key
        # ("" is not a valid identifier and was almost certainly a typo)
        assert "total mers" in c.data["correction"]
        assert "" not in c.data["correction"]

    def test_set_overlap_filtering(self, canu_dir):
        c = CanuScanner(canu_dir)
        c.set_overlap_filtering()
        assert "overlap filtering stats text" in c.data["correction"]["overlap filtering"]

    def test_set_read_correction(self, canu_dir):
        c = CanuScanner(canu_dir)
        c.set_read_correction()
        assert "correction summary text" in c.data["correction"]["read correction"]

    def test_plot_correction_check1_handles_missing_files_gracefully(self, canu_dir):
        """tn/fn/fp/tp log files are absent from the fixture; the method must
        not raise (real Canu runs often lack some of these), while still
        raising on genuinely unexpected errors."""
        c = CanuScanner(canu_dir)
        c.plot_correction_check1()  # should not raise

    def test_hist_read_length2(self, canu_dir):
        c = CanuScanner(canu_dir)
        df = c.hist_read_length2()
        assert len(df) == 2
