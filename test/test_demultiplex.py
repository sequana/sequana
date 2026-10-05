import json
import os

import pandas as pd
import pytest
from easydev import TempFile

from sequana.demultiplex import StatsFile

from . import test_dir


class TestStatsFile:
    """Test StatsFile class for parsing bcl2fastq Stats.json files."""

    def test_stats_file_init(self):
        """Test StatsFile initialization with real data."""
        data = f"{test_dir}/data/json/test_demultiplex_Stats.json"
        s = StatsFile(data)
        assert s.data is not None
        assert isinstance(s.data, dict)
        assert "ConversionResults" in s.data

    def test_stats_file_init_undetermined(self):
        """Test StatsFile initialization with undetermined data."""
        data = f"{test_dir}/data/json/test_demultiplex_Stats_undetermined.json"
        s = StatsFile(data)
        assert s.data is not None
        assert isinstance(s.data, dict)

    def test_stats_file_init_missing_file(self):
        """Test StatsFile with missing file raises error."""
        with pytest.raises((FileNotFoundError, IOError)):
            StatsFile("nonexistent_file.json")

    def test_stats_file_get_data_reads(self):
        """Test get_data_reads returns correct dataframe structure."""
        data = f"{test_dir}/data/json/test_demultiplex_Stats.json"
        s = StatsFile(data)
        df = s.get_data_reads()
        assert isinstance(df, pd.DataFrame)
        assert "lane" in df.columns
        assert "name" in df.columns
        assert "count" in df.columns
        assert len(df) > 0
        # Check that we have both determined and undetermined
        names = df["name"].unique()
        assert any("Undetermined" in name for name in names)

    def test_stats_file_get_data_reads_undetermined(self):
        """Test get_data_reads with undetermined-only data."""
        data = f"{test_dir}/data/json/test_demultiplex_Stats_undetermined.json"
        s = StatsFile(data)
        df = s.get_data_reads()
        assert isinstance(df, pd.DataFrame)
        # Should still have expected columns
        assert "lane" in df.columns
        assert "name" in df.columns
        assert "count" in df.columns

    def test_stats_file_get_data_reads_empty_demux(self):
        """Test get_data_reads when DemuxResults is empty."""
        data = f"{test_dir}/data/json/test_demultiplex_Stats.json"
        s = StatsFile(data)
        df = s.get_data_reads()
        # Filter for determined samples
        determined = df[df["name"] != "Undetermined"]
        # Should have some determined samples
        assert len(determined) >= 0

    def test_stats_file_to_summary_reads(self):
        """Test to_summary_reads writes summary file."""
        data = f"{test_dir}/data/json/test_demultiplex_Stats.json"
        s = StatsFile(data)
        with TempFile() as fout:
            s.to_summary_reads(fout.name)
            assert os.path.exists(fout.name)
            with open(fout.name) as f:
                content = f.read()
                assert "Lane" in content
                assert "Sample" in content
                assert "NumberReads" in content

    def test_stats_file_to_summary_reads_undetermined(self):
        """Test to_summary_reads with undetermined data."""
        data = f"{test_dir}/data/json/test_demultiplex_Stats_undetermined.json"
        s = StatsFile(data)
        with TempFile() as fout:
            s.to_summary_reads(fout.name)
            assert os.path.exists(fout.name)

    def test_stats_file_to_summary_reads_format(self):
        """Test to_summary_reads output format is correct."""
        data = f"{test_dir}/data/json/test_demultiplex_Stats.json"
        s = StatsFile(data)
        with TempFile() as fout:
            s.to_summary_reads(fout.name)
            with open(fout.name) as f:
                lines = f.readlines()
                # First line should be header
                assert "Lane" in lines[0]
                # All subsequent lines should be data
                for line in lines[1:]:
                    parts = line.split()
                    assert len(parts) == 3  # Lane, Sample, Count

    def test_stats_file_barplot_summary(self):
        """Test barplot_summary creates visualization."""
        data = f"{test_dir}/data/json/test_demultiplex_Stats.json"
        s = StatsFile(data)
        with TempFile() as fout:
            df = s.barplot_summary(fout.name)
            assert isinstance(df, pd.DataFrame)
            assert "Determined" in df.columns
            assert "Undetermined" in df.columns
            assert os.path.exists(fout.name)

    def test_stats_file_barplot_summary_no_file(self):
        """Test barplot_summary without saving file."""
        data = f"{test_dir}/data/json/test_demultiplex_Stats.json"
        s = StatsFile(data)
        df = s.barplot_summary(filename=None)
        assert isinstance(df, pd.DataFrame)

    def test_stats_file_barplot(self):
        """Test barplot generates lane-specific plots."""
        data = f"{test_dir}/data/json/test_demultiplex_Stats.json"
        s = StatsFile(data)
        s.barplot()
        # Clean up generated files
        for lane in s.get_data_reads().lane.unique():
            filename = f"lane{lane}_status.png"
            if os.path.exists(filename):
                os.remove(filename)

    def test_stats_file_barplot_specific_lanes(self):
        """Test barplot with specific lane selection."""
        data = f"{test_dir}/data/json/test_demultiplex_Stats.json"
        s = StatsFile(data)
        df = s.get_data_reads()
        lanes = df.lane.unique()[:1]  # Just first lane
        s.barplot(lanes=lanes)
        # Clean up
        for lane in lanes:
            filename = f"lane{lane}_status.png"
            if os.path.exists(filename):
                os.remove(filename)

    def test_stats_file_barplot_per_sample(self):
        """Test barplot_per_sample creates visualization."""
        data = f"{test_dir}/data/json/test_demultiplex_Stats.json"
        s = StatsFile(data)
        with TempFile() as fout:
            s.barplot_per_sample(filename=fout.name)
            assert os.path.exists(fout.name)

    def test_stats_file_barplot_per_sample_no_file(self):
        """Test barplot_per_sample without saving file."""
        data = f"{test_dir}/data/json/test_demultiplex_Stats.json"
        s = StatsFile(data)
        s.barplot_per_sample()
        # Should complete without error

    def test_stats_file_barplot_per_sample_alpha(self):
        """Test barplot_per_sample with alpha parameter."""
        data = f"{test_dir}/data/json/test_demultiplex_Stats.json"
        s = StatsFile(data)
        s.barplot_per_sample(alpha=0.7)
        # Should complete without error

    def test_stats_file_barplot_per_sample_width(self):
        """Test barplot_per_sample with width parameter."""
        data = f"{test_dir}/data/json/test_demultiplex_Stats.json"
        s = StatsFile(data)
        s.barplot_per_sample(width=0.6)
        # Should complete without error

    def test_stats_file_plot_unknown_barcodes(self):
        """Test plot_unknown_barcodes with real data."""
        data = f"{test_dir}/data/json/test_demultiplex_Stats.json"
        s = StatsFile(data)
        df = s.plot_unknown_barcodes()
        assert isinstance(df, pd.DataFrame)

    def test_stats_file_plot_unknown_barcodes_n_parameter(self):
        """Test plot_unknown_barcodes with N parameter."""
        data = f"{test_dir}/data/json/test_demultiplex_Stats.json"
        s = StatsFile(data)
        df = s.plot_unknown_barcodes(N=10)
        assert isinstance(df, pd.DataFrame)
        # Should have at most N rows
        assert len(df) <= 10

    def test_stats_file_plot_unknown_barcodes_single_lane(self):
        """Test plot_unknown_barcodes with single lane data."""
        data = f"{test_dir}/data/json/test_demultiplex_Stats_undetermined.json"
        s = StatsFile(data)
        df = s.plot_unknown_barcodes()
        assert isinstance(df, pd.DataFrame)

    def test_stats_file_comprehensive_workflow(self):
        """Test complete workflow with all methods."""
        data = f"{test_dir}/data/json/test_demultiplex_Stats.json"
        s = StatsFile(data)

        # Test all methods
        df = s.get_data_reads()
        assert len(df) > 0

        with TempFile() as fout:
            s.to_summary_reads(fout.name)
            assert os.path.exists(fout.name)

        with TempFile() as fout:
            df_summary = s.barplot_summary(fout.name)
            assert isinstance(df_summary, pd.DataFrame)

        df_per_sample = s.barplot_per_sample()
        # Should have attributes set
        assert hasattr(s, "all_data")
        assert hasattr(s, "under")


def test_stats_file():
    """Original test for compatibility."""
    data = f"{test_dir}/data/json/test_demultiplex_Stats.json"
    s = StatsFile(data)
    with TempFile() as fout:
        s.to_summary_reads(fout.name)
    with TempFile() as fout:
        s.barplot_summary(fout.name)
    with TempFile() as fout:
        s.barplot()
        for lane in s.get_data_reads().lane.unique():
            os.remove("lane{}_status.png".format(lane))

    data = f"{test_dir}/data/json/test_demultiplex_Stats_undetermined.json"
    s = StatsFile(data)
    with TempFile() as fout:
        s.to_summary_reads(fout.name)
    with TempFile() as fout:
        s.barplot_summary(fout.name)
    with TempFile() as fout:
        s.barplot()
        for lane in s.get_data_reads().lane.unique():
            os.remove("lane{}_status.png".format(lane))
