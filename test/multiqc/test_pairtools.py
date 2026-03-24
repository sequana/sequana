import io
import math

import pytest

try:
    from sequana.multiqc.pairtools import (
        cumsums_to_rangesums,
        edges_to_intervals,
        geometric_mean,
        humanize_genomic_dist,
        humanize_genomic_interval,
        parse_dist_range,
        read_stats_from_file,
    )

    # --- humanize_genomic_dist ---

    def test_humanize_genomic_dist_bp():
        assert humanize_genomic_dist(500) == "500bp"

    def test_humanize_genomic_dist_kb():
        assert humanize_genomic_dist(5000) == "5Kb"

    def test_humanize_genomic_dist_mb():
        assert humanize_genomic_dist(3_000_000) == "3Mb"

    def test_humanize_genomic_dist_with_units():
        assert humanize_genomic_dist(2, units=1000) == "2Kb"

    # --- humanize_genomic_interval ---

    def test_humanize_genomic_interval_last_open():
        assert humanize_genomic_interval(5_000_000, None) == ">5Mb"

    def test_humanize_genomic_interval_first():
        assert humanize_genomic_interval(0, 1000) == "<1Kb"

    def test_humanize_genomic_interval_both_mb():
        assert humanize_genomic_interval(2_000_000, 4_000_000) == "2-4Mb"

    def test_humanize_genomic_interval_both_kb():
        assert humanize_genomic_interval(2_000, 5_000) == "2-5Kb"

    def test_humanize_genomic_interval_end_kb():
        assert humanize_genomic_interval(500, 2_000) == "0.5-2Kb"

    def test_humanize_genomic_interval_bp():
        assert humanize_genomic_interval(100, 900) == "100-900bp"

    def test_humanize_genomic_interval_end_mb():
        result = humanize_genomic_interval(500_000, 2_000_000)
        assert "Mb" in result

    # --- edges_to_intervals ---

    def test_edges_to_intervals_without_zero():
        assert edges_to_intervals([1, 2, 3]) == [(0, 1), (1, 2), (2, 3), (3, None)]

    def test_edges_to_intervals_with_zero():
        assert edges_to_intervals([0, 1, 2, 3]) == [(0, 1), (1, 2), (2, 3), (3, None)]

    def test_edges_to_intervals_single():
        assert edges_to_intervals([5]) == [(0, 5), (5, None)]

    # --- cumsums_to_rangesums ---

    def test_cumsums_to_rangesums_basic():
        assert cumsums_to_rangesums([80, 50, 10], total=100) == [20, 30, 40, 10]

    def test_cumsums_to_rangesums_single():
        assert cumsums_to_rangesums([60], total=100) == [40, 60]

    # --- parse_dist_range ---

    def test_parse_dist_range_open_ended():
        assert parse_dist_range("10000+") == (10000, None)

    def test_parse_dist_range_interval():
        assert parse_dist_range("100-2000") == (100, 2000)

    def test_parse_dist_range_zero_start():
        assert parse_dist_range("0-10") == (0, 10)

    def test_parse_dist_range_invalid():
        with pytest.raises(ValueError, match="cannot parse"):
            parse_dist_range("invalid")

    # --- geometric_mean ---

    def test_geometric_mean_zero_start():
        assert geometric_mean(0, 100) == 0

    def test_geometric_mean_open_ended():
        assert geometric_mean(500, None) == 500

    def test_geometric_mean_normal():
        assert abs(geometric_mean(4, 9) - 6.0) < 1e-10

    def test_geometric_mean_equal():
        assert abs(geometric_mean(16, 16) - 16.0) < 1e-10

    # --- read_stats_from_file ---

    MINIMAL_STATS = """\
total\t1000
total_unmapped\t100
total_single_sided_mapped\t50
total_mapped\t850
total_dups\t200
total_nodups\t650
cis\t400
trans\t250
pair_types/UU\t600
pair_types/NU\t100
"""

    def test_read_stats_from_file_basic():
        result = read_stats_from_file(io.StringIO(MINIMAL_STATS))
        assert result is not None
        assert result["total"] == 1000
        assert result["total_unmapped"] == 100
        assert abs(result["frac_unmapped"] - 10.0) < 1e-6
        assert abs(result["frac_mapped"] - 85.0) < 1e-6
        assert abs(result["frac_cis"] - 400 / 650 * 100) < 1e-6
        assert result["pair_types"]["UU"] == 600

    def test_read_stats_from_file_missing_required_key():
        incomplete = """\
total\t1000
total_unmapped\t100
total_single_sided_mapped\t50
total_mapped\t850
total_dups\t200
total_nodups\t650
trans\t250
"""
        result = read_stats_from_file(io.StringIO(incomplete))
        assert result is None

    def test_read_stats_from_file_with_cis_dist():
        stats_with_cis = MINIMAL_STATS + "cis_1kb+\t300\ncis_10kb+\t150\n"
        result = read_stats_from_file(io.StringIO(stats_with_cis))
        assert result is not None
        assert result["cis_dist"] is not None
        assert "trans" in result["cis_dist"]

    def test_read_stats_from_file_skips_non_numeric():
        stats_with_text = MINIMAL_STATS + "some_key\tnot_a_number\n"
        result = read_stats_from_file(io.StringIO(stats_with_text))
        assert result is not None

    def test_read_stats_from_file_with_dist_freq():
        dist_freq_lines = ""
        for orient in ["++", "+-", "-+", "--"]:
            dist_freq_lines += f"dist_freq/100-500/{orient}\t50\n"
            dist_freq_lines += f"dist_freq/500-1000/{orient}\t30\n"
        result = read_stats_from_file(io.StringIO(MINIMAL_STATS + dist_freq_lines))
        assert result is not None
        assert result["dist_freq"] is not None
        assert "all" in result["dist_freq"]

except ImportError:
    pass
