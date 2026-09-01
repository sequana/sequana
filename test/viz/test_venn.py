import pytest

from sequana.viz.venn import get_layout, plot_venn


def test3():
    A = set([1, 2, 3, 4, 5, 6, 7, 8, 9])
    B = set([6, 7, 8, 9, 10, 11, 12, 13])
    C = set([4, 5, 6, 7, 8, 9])
    plot_venn((A, B, C), labels=("A", "B", "C"))
    plot_venn((A, B, C), labels=("A", "B", "C"), weighted=True)


def test2():
    A = set([1, 2, 3, 4, 5, 6, 7, 8, 9])
    B = set([6, 7, 8, 9, 10, 11, 12, 13])
    plot_venn((A, B), labels=("A", "B"))
    plot_venn((A, B), labels=("A", "B"), weighted=True)


@pytest.mark.parametrize("num_sets", [2, 3, 4, 5, 6])
def test_get_layout(num_sets):
    layout = get_layout(num_sets)
    assert len(layout["shapes"]) == num_sets, f"Expected {num_sets} shapes for {num_sets}-set layout"
    assert len(layout["names"]) == num_sets, f"Expected {num_sets} name positions for {num_sets}-set layout"
    for shape_type, params in layout["shapes"]:
        assert shape_type == "ellipse"
        assert "xy" in params and "width" in params and "height" in params and "angle" in params
