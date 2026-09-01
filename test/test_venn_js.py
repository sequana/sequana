import json
import re

import pandas as pd

from sequana.utils.venn_js import DynamicVenn


def test_encode_log2_fold_change():
    l2fc = pd.Series([0.1, 1.0, -1.0, 3.2, -8, None])
    significant = pd.Series([True, True, True, False, True, True])

    codes = DynamicVenn.encode_log2_fold_change(l2fc, significant, step=0.5, nbins=9)

    # 0.1 is in the first bin, 1.0 in the third one, and the sign gives the
    # direction of the regulation
    assert list(codes) == [1, 3, -3, 0, -9, 0]

    # the last bin gathers everything above (nbins - 1) * step
    assert abs(codes[4]) == 9


def test_dynamic_venn_html():
    names = ["A_vs_B", "C_vs_D", "E_vs_F"]
    rows = [[1, 0, 0, 10], [1, 2, 0, 5], [0, 0, -3, 2]]

    html = DynamicVenn(names, rows, step=0.25, nbins=17, html_id="myvenn").to_html()

    # the controls and the container the javascript writes into
    for name in names:
        assert name in html
    assert 'id="myvenn_plot"' in html
    assert 'id="myvenn_slider"' in html
    assert html.count('type="checkbox"') == len(names)

    # all/up/down and the mode where a comparison gives an up and a down set
    assert html.count('type="radio"') == 4
    assert 'value="split"' in html

    # the data itself, including the layouts of the 2 and 3-set diagrams
    data = json.loads(re.search(r"var D = (\{.*?\});\n", html, re.DOTALL).group(1))
    assert data["rows"] == rows
    assert data["names"] == names
    assert sorted(data["layouts"]) == ["2", "3"]
    assert data["selected"] == [0, 1, 2]


def test_dynamic_venn_max_sets():
    names = [f"c{i}" for i in range(8)]
    html = DynamicVenn(names, [[1] * 8 + [3]], html_id="v").to_html()

    data = json.loads(re.search(r"var D = (\{.*?\});\n", html, re.DOTALL).group(1))

    # only 6 sets can be drawn, so only 6 comparisons are selected by default
    assert data["selected"] == list(range(6))
    assert sorted(data["layouts"]) == ["2", "3", "4", "5", "6"]

    # layouts for 4-6 sets must have non-empty shapes and names
    for n in range(2, 7):
        layout = data["layouts"][str(n)]
        assert len(layout["shapes"]) == n, f"Layout for {n} sets should have {n} shapes"
        assert len(layout["names"]) == n, f"Layout for {n} sets should have {n} name positions"
