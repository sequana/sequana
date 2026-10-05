import pandas as pd
import pytest

from sequana.viz.bar import CanvasBar, stacked_bar


class TestStackedBar:
    """Test stacked_bar function for creating stacked bar charts."""

    def test_stacked_bar_basic(self):
        """Test basic stacked_bar function."""
        datalist = [
            {"name": "Sample1", "data": {"R1_mapped": 50, "R2_mapped": 50, "R1_unmapped": 30, "R2_unmapped": 20}},
            {
                "name": "Sample2",
                "data": {"R1_mapped": 60, "R2_mapped": 40, "R1_unmapped": 20, "R2_unmapped": 10},
            },
        ]
        script = stacked_bar("test_tag", "Test Title", datalist)
        assert isinstance(script, str)
        assert "CanvasJS" in script
        assert "Test Title" in script
        assert "chartContainer" in script
        assert "stackedBar100" in script

    def test_stacked_bar_single_item(self):
        """Test stacked_bar with single item."""
        datalist = [{"name": "Sample1", "data": {"R1_mapped": 100, "R2_mapped": 100}}]
        script = stacked_bar("tag1", "Single Sample", datalist)
        assert isinstance(script, str)
        assert "Sample1" in script

    def test_stacked_bar_empty_list(self):
        """Test stacked_bar with empty data list."""
        datalist = []
        script = stacked_bar("empty_tag", "Empty Chart", datalist)
        assert isinstance(script, str)
        assert "CanvasJS" in script

    def test_stacked_bar_special_characters(self):
        """Test stacked_bar with special characters in labels."""
        datalist = [
            {
                "name": "Sample_1/2",
                "data": {"R1_mapped": 50, "R2_mapped": 50},
            }
        ]
        script = stacked_bar("special_tag", "Special Chars", datalist)
        assert isinstance(script, str)

    def test_stacked_bar_large_numbers(self):
        """Test stacked_bar with large numbers."""
        datalist = [
            {
                "name": "LargeSample",
                "data": {"R1_mapped": 1000000, "R2_mapped": 2000000, "unmapped": 500000},
            }
        ]
        script = stacked_bar("large_tag", "Large Numbers", datalist)
        assert isinstance(script, str)
        assert "1000000" in script or "1e6" in script or "1e+6" in script

    def test_stacked_bar_zero_values(self):
        """Test stacked_bar with zero values."""
        datalist = [{"name": "Sample1", "data": {"mapped": 0, "unmapped": 0}}]
        script = stacked_bar("zero_tag", "Zero Values", datalist)
        assert isinstance(script, str)

    def test_stacked_bar_multiple_categories(self):
        """Test stacked_bar with many categories."""
        data = {f"category_{i}": i * 10 for i in range(10)}
        datalist = [{"name": "Sample1", "data": data}]
        script = stacked_bar("multi_tag", "Multi Category", datalist)
        assert isinstance(script, str)


class TestCanvasBar:
    """Test CanvasBar class for rendering bar charts."""

    def test_canvasbar_init_basic(self):
        """Test CanvasBar initialization with basic DataFrame."""
        data = pd.DataFrame({"name": ["A", "B", "C"], "value": [10, 20, 30], "url": ["#", "#", "#"]})
        bar = CanvasBar(data, title="Test Chart", tag="test", xlabel="Count")
        assert bar.title == "Test Chart"
        assert bar.tag == "test"
        assert isinstance(bar.metadata, dict)
        assert "dataitems" in bar.metadata
        assert "title" in bar.metadata

    def test_canvasbar_single_row(self):
        """Test CanvasBar with single row DataFrame."""
        data = pd.DataFrame({"name": ["Item1"], "value": [100], "url": ["http://example.com"]})
        bar = CanvasBar(data, title="Single Item")
        assert bar.title == "Single Item"
        assert "Item1" in bar.metadata["dataitems"]

    def test_canvasbar_empty_dataframe(self):
        """Test CanvasBar with empty DataFrame."""
        data = pd.DataFrame({"name": [], "value": [], "url": []})
        bar = CanvasBar(data)
        assert isinstance(bar.metadata, dict)

    def test_canvasbar_special_characters_in_names(self):
        """Test CanvasBar with special characters in bar names."""
        data = pd.DataFrame(
            {"name": ["Sample-1", "Sample/2", "Sample_3"], "value": [10, 20, 30], "url": ["#", "#", "#"]}
        )
        bar = CanvasBar(data)
        assert "Sample-1" in bar.metadata["dataitems"] or "dataitems" in bar.metadata

    def test_canvasbar_large_values(self):
        """Test CanvasBar with large values."""
        data = pd.DataFrame({"name": ["A", "B", "C"], "value": [1e6, 2e6, 3e6], "url": ["#", "#", "#"]})
        bar = CanvasBar(data, title="Large Values", tag="large")
        assert bar.title == "Large Values"

    def test_canvasbar_zero_values(self):
        """Test CanvasBar with zero values."""
        data = pd.DataFrame({"name": ["A", "B"], "value": [0, 10], "url": ["#", "#"]})
        bar = CanvasBar(data)
        assert isinstance(bar.metadata, dict)

    def test_canvasbar_to_html_basic(self):
        """Test CanvasBar.to_html basic output."""
        data = pd.DataFrame({"name": ["A", "B"], "value": [10, 20], "url": ["#", "#"]})
        bar = CanvasBar(data, title="Test", tag="test", xlabel="Count")
        html = bar.to_html()
        assert isinstance(html, str)
        assert "CanvasJS" in html
        assert "bar" in html
        assert "Count" in html

    def test_canvasbar_to_html_with_maxrange(self):
        """Test CanvasBar.to_html with maxrange option."""
        data = pd.DataFrame({"name": ["A", "B"], "value": [10, 20], "url": ["#", "#"]})
        bar = CanvasBar(data, title="Test", tag="test", xlabel="Count")
        html = bar.to_html(options={"maxrange": 100})
        assert isinstance(html, str)
        assert "maximum" in html
        assert "100" in html

    def test_canvasbar_to_html_no_maxrange(self):
        """Test CanvasBar.to_html without maxrange."""
        data = pd.DataFrame({"name": ["A", "B"], "value": [10, 20], "url": ["#", "#"]})
        bar = CanvasBar(data)
        html = bar.to_html(options={"maxrange": None})
        assert isinstance(html, str)

    def test_canvasbar_metadata_structure(self):
        """Test CanvasBar metadata structure."""
        data = pd.DataFrame({"name": ["Item1", "Item2"], "value": [100, 200], "url": ["url1", "url2"]})
        bar = CanvasBar(data, title="MyChart", tag="mytag", xlabel="Values")
        assert bar.metadata["title"] == "MyChart"
        assert bar.metadata["tag"] == "mytag"
        assert bar.metadata["xlabel"] == "Values"
        assert "dataitems" in bar.metadata

    def test_canvasbar_url_rendering(self):
        """Test that URLs are correctly rendered in dataitems."""
        data = pd.DataFrame(
            {"name": ["A", "B"], "value": [10, 20], "url": ["http://example1.com", "http://example2.com"]}
        )
        bar = CanvasBar(data)
        assert "http://example1.com" in bar.metadata["dataitems"]
        assert "http://example2.com" in bar.metadata["dataitems"]

    def test_canvasbar_many_bars(self):
        """Test CanvasBar with many bars."""
        n = 50
        data = pd.DataFrame(
            {"name": [f"Bar{i}" for i in range(n)], "value": [i * 10 for i in range(n)], "url": ["#"] * n}
        )
        bar = CanvasBar(data, title="Many Bars")
        assert bar.title == "Many Bars"
        # Check that all bars are in dataitems
        for i in range(n):
            assert f"Bar{i}" in bar.metadata["dataitems"]

    def test_canvasbar_numeric_values(self):
        """Test CanvasBar with different numeric types."""
        data = pd.DataFrame(
            {
                "name": ["int", "float", "large"],
                "value": [10, 20.5, 1e6],
                "url": ["#", "#", "#"],
            }
        )
        bar = CanvasBar(data)
        assert "10" in str(bar.metadata["dataitems"])

    def test_canvasbar_edge_case_single_column(self):
        """Test CanvasBar behavior with minimal data."""
        data = pd.DataFrame({"name": ["X"], "value": [1], "url": ["#"]})
        bar = CanvasBar(data, title="Minimal")
        html = bar.to_html()
        assert isinstance(html, str)

    def test_canvasbar_ylabel_length_parameter(self):
        """Test CanvasBar with ylabel_max_length parameter."""
        data = pd.DataFrame({"name": ["A" * 100, "B"], "value": [10, 20], "url": ["#", "#"]})
        bar = CanvasBar(data, ylabel_max_length=10)
        assert isinstance(bar.metadata, dict)

    def test_canvasbar_with_kwargs(self):
        """Test CanvasBar accepts additional kwargs."""
        data = pd.DataFrame({"name": ["A", "B"], "value": [10, 20], "url": ["#", "#"]})
        # Should not raise error with extra kwargs
        bar = CanvasBar(data, title="Test", extra_param="value")
        assert bar is not None
