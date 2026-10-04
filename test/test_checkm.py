"""Tests for sequana.checkm."""
import pytest

from sequana.checkm import CheckM, MultiCheckM

# CheckM's lineage_wf report format: a fixed-width table using >=2 spaces as
# field separator; the record row is always line index 3 (0-indexed) after
# 3 header/separator lines.
CHECKM_REPORT = (
    "--------------------------------------------------------------------------------\n"
    "  Bin Id      Marker lineage      # genomes   # markers   # marker sets    0    1    2   3   4   5+   Completeness   Contamination   Strain heterogeneity\n"
    "--------------------------------------------------------------------------------\n"
    "  sample1  k__Bacteria (UID203)      5449         104           58         0   103   1   0   0   0       99.14           1.72               0.00\n"
    "--------------------------------------------------------------------------------\n"
)


@pytest.fixture
def checkm_report_file(tmpdir):
    path = tmpdir.join("results.txt")
    path.write(CHECKM_REPORT)
    return str(path)


class TestCheckM:
    def test_parses_sample_name(self, checkm_report_file):
        c = CheckM(checkm_report_file)
        assert c.df["sample"] == "sample1"

    def test_parses_marker_lineage(self, checkm_report_file):
        c = CheckM(checkm_report_file)
        assert c.df["marker_lineage"] == "k__Bacteria (UID203)"

    def test_parses_numeric_fields(self, checkm_report_file):
        c = CheckM(checkm_report_file)
        assert c.df["Completeness"] == 99.14
        assert c.df["Contamination"] == 1.72
        assert c.df["Strain heterogeneity"] == 0.0
        assert c.df["#genomes"] == 5449.0

    def test_index_matches_header(self, checkm_report_file):
        c = CheckM(checkm_report_file)
        expected = [
            "sample",
            "marker_lineage",
            "#genomes",
            "#markers",
            "#marker_sets",
            "0",
            "1",
            "2",
            "3",
            "4",
            "5+",
            "Completeness",
            "Contamination",
            "Strain heterogeneity",
        ]
        assert list(c.df.index) == expected

    def test_missing_file_raises(self):
        with pytest.raises(FileNotFoundError):
            CheckM("/tmp/this_file_should_not_exist_checkm.txt")


class TestMultiCheckM:
    def test_concatenates_multiple_samples(self, tmpdir, checkm_report_file):
        second = tmpdir.join("results2.txt")
        second.write(CHECKM_REPORT)

        m = MultiCheckM([checkm_report_file, str(second)])
        assert m.df.shape[1] == 2
        assert (m.df.loc["sample"] == "sample1").all()

    def test_skips_unreadable_files_with_warning(self, checkm_report_file, caplog):
        m = MultiCheckM([checkm_report_file, "/tmp/this_file_should_not_exist_checkm.txt"])
        # only the valid file contributes a column
        assert m.df.shape[1] == 1

    def test_all_files_missing_gives_empty_concat(self):
        m = MultiCheckM(["/tmp/missing_a.txt", "/tmp/missing_b.txt"])
        assert m.df.empty
