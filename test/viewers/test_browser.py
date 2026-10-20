"""Interaction tests of the generated page in a real browser (needs playwright + chromium)."""
import re

import numpy as np
import pytest

from sequana.viewers.genome_track import build_page, write_page

sync_api = pytest.importorskip("playwright.sync_api")
pytestmark = pytest.mark.slow

LENGTH = 300_000


def synthetic_page(tmp_path_factory, with_gc=True, with_annotation=True, with_thresholds=True):
    rng = np.random.default_rng(0)
    signal = rng.poisson(40, LENGTH).astype(float)
    signal[100_000:100_400] = 5  # a deletion
    baseline = np.full(LENGTH, 40.0)
    regions = [
        {"start": 100_000, "end": 100_400, "mean_signal": 5.0, "mean_baseline": 40.0, "score": -9.0, "type": "low"},
        {"start": 200_000, "end": 200_003, "mean_signal": 90.0, "mean_baseline": 40.0, "score": 7.0, "type": "high"},
    ]
    features = [
        {
            "start": 99_000,
            "end": 101_000,
            "strand": "+",
            "type": "gene",
            "name": "geneA",
            "product": "a kinase",
            "id": "A",
            "lane": 0,
        },
        {
            "start": 150_000,
            "end": 153_000,
            "strand": "-",
            "type": "gene",
            "name": "geneB",
            "product": "a transporter",
            "id": "B",
            "lane": 0,
        },
    ]
    kwargs = {}
    if with_gc:
        kwargs["gc"] = np.clip(0.5 + 0.1 * np.sin(np.arange(LENGTH) / 3000), 0, 1)
    if with_annotation:
        kwargs.update(features=features, lanes=1)
    if with_thresholds:
        kwargs["thresholds"] = {"low": -4, "high": 4, "mu": 1.0, "sigma": 0.15}
    path = tmp_path_factory.mktemp("viewer") / "page.html"
    write_page(path, "synthetic", signal, baseline, regions, signal_mean=40, **kwargs)
    return path


@pytest.fixture(scope="module")
def browser():
    with sync_api.sync_playwright() as p:
        try:
            b = p.chromium.launch()
        except Exception as err:  # browser not installed
            pytest.skip(f"chromium not available: {err}")
        yield b
        b.close()


@pytest.fixture
def open_page(browser, tmp_path_factory):
    pages = []

    def _open(**kwargs):
        path = synthetic_page(tmp_path_factory, **kwargs)
        page = browser.new_page(viewport={"width": 1500, "height": 1100})
        page.errors = []
        page.on("pageerror", lambda e: page.errors.append(str(e)))
        page.on("console", lambda m: page.errors.append(m.text) if m.type == "error" else None)
        page.goto(path.as_uri())
        page.wait_for_function("window.__ready === true", timeout=30000)
        page.wait_for_timeout(200)
        pages.append(page)
        return page

    yield _open
    for page in pages:
        page.close()


def view(page):
    m = re.search(r"view ([\d,]+) - ([\d,]+) bp", page.inner_text("#status"))
    return int(m.group(1).replace(",", "")), int(m.group(2).replace(",", ""))


def yrange(page):
    return page.evaluate("() => [window.__viewer.main.scales.y.min, window.__viewer.main.scales.y.max]")


def test_loads_without_errors_and_shows_everything(open_page):
    page = open_page()
    assert view(page) == (1, LENGTH)
    assert not page.errors
    assert page.locator("#roitable tbody tr").count() == 2
    assert page.locator("#annplot canvas").count() == 1


def test_wheel_zoom_keeps_pointer_position_and_refits_y(open_page):
    page = open_page()
    box = page.locator("#plot .u-over").bounding_box()
    x = box["x"] + box["width"] * 0.4
    page.mouse.move(x, box["y"] + 100)
    before = page.evaluate(
        "(x) => window.__viewer.main.posToVal(x - document.querySelector('#plot .u-over').getBoundingClientRect().left, 'x')",
        x,
    )
    for _ in range(4):
        page.mouse.wheel(0, -120)
        page.wait_for_timeout(100)
    page.wait_for_timeout(300)
    lo, hi = view(page)
    assert hi - lo < 0.4 * LENGTH
    after = page.evaluate(
        "(x) => window.__viewer.main.posToVal(x - document.querySelector('#plot .u-over').getBoundingClientRect().left, 'x')",
        x,
    )
    assert abs(after - before) < 0.02 * (hi - lo) + 50
    drawn = page.evaluate("() => { const d = window.__viewer.main.data[0]; return [d[0], d[d.length - 1]]; }")
    assert drawn[0] <= lo and drawn[1] >= hi  # data cover the whole window
    for _ in range(30):
        page.mouse.wheel(0, 120)
        page.wait_for_timeout(30)
    page.wait_for_timeout(300)
    assert view(page) == (1, LENGTH)
    assert not page.errors


def test_rectangle_zoom_sets_x_and_y(open_page):
    page = open_page()
    box = page.locator("#plot .u-over").bounding_box()
    y0 = yrange(page)
    xa, xb = box["x"] + box["width"] * 0.3, box["x"] + box["width"] * 0.4
    ya, yb = box["y"] + box["height"] * 0.5, box["y"] + box["height"] * 0.8
    page.mouse.move(xa, ya)
    page.mouse.down()
    page.mouse.move(xb, yb, steps=6)
    page.mouse.up()
    page.wait_for_timeout(400)
    lo, hi = view(page)
    assert hi - lo < 0.15 * LENGTH
    top = y0[1] - (ya - box["y"]) / box["height"] * (y0[1] - y0[0])
    assert yrange(page)[1] == pytest.approx(top, abs=0.06 * (y0[1] - y0[0]))
    assert not page.errors


def test_right_button_and_shift_drag_pan(open_page):
    page = open_page()
    box = page.locator("#plot .u-over").bounding_box()
    for _ in range(3):
        page.mouse.move(box["x"] + 500, box["y"] + 100)
        page.mouse.wheel(0, -120)
        page.wait_for_timeout(100)
    page.wait_for_timeout(300)
    start = view(page)
    page.mouse.move(box["x"] + 600, box["y"] + 100)
    page.mouse.down(button="right")
    page.mouse.move(box["x"] + 300, box["y"] + 100, steps=6)
    page.mouse.up(button="right")
    page.wait_for_timeout(400)
    moved = view(page)
    assert moved[0] > start[0] and abs((moved[1] - moved[0]) - (start[1] - start[0])) <= 2  # same width
    page.keyboard.down("Shift")
    page.mouse.move(box["x"] + 300, box["y"] + 100)
    page.mouse.down()
    page.mouse.move(box["x"] + 500, box["y"] + 100, steps=6)
    page.mouse.up()
    page.keyboard.up("Shift")
    page.wait_for_timeout(400)
    assert view(page)[0] < moved[0]
    assert not page.errors


def test_clicking_a_region_zooms_to_it_and_gene_names_are_listed(open_page):
    page = open_page()
    assert page.locator("#roitable th", has_text="Gene").count() == 1
    first = page.locator("#roitable tbody tr").first  # sorted by size: the deletion
    assert "geneA" in first.inner_text() and "a kinase" in first.inner_text()
    assert "intergenic" in page.locator("#roitable tbody tr").nth(1).inner_text()
    first.click()
    page.wait_for_timeout(600)
    lo, hi = view(page)
    assert lo <= 100_000 and hi >= 100_400 and hi - lo < 10_000
    page.fill("#search", "transporter")
    page.wait_for_timeout(200)
    assert page.locator("#roitable tbody tr").count() == 0  # geneB has no region
    page.fill("#search", "kinase")
    page.wait_for_timeout(200)
    assert page.locator("#roitable tbody tr").count() == 1


def test_gene_track_hover_and_click(open_page):
    page = open_page()
    page.evaluate("() => window.__viewer.main.setScale('x', {min: 140000, max: 165000})")
    page.wait_for_timeout(500)
    box = page.locator("#annplot .u-over").bounding_box()
    x = box["x"] + (151_500 - 140_000) / 25_000 * box["width"]
    y = box["y"] + 8
    page.mouse.move(x, y)
    page.wait_for_timeout(200)
    assert page.is_visible("#tip") and "geneB" in page.inner_text("#tip")
    page.mouse.click(x, y)
    page.wait_for_timeout(500)
    lo, hi = view(page)
    assert lo <= 150_000 and hi >= 153_000 and hi - lo < 10_000


def test_gc_axis_is_fixed_and_aligned(open_page):
    page = open_page()
    page.evaluate("() => window.__viewer.main.setScale('x', {min: 50000, max: 70000})")
    page.wait_for_timeout(400)
    assert page.evaluate("() => [window.__viewer.main.scales.gc.min, window.__viewer.main.scales.gc.max]") == [0, 1]
    page.check("#gcpanel")
    page.wait_for_timeout(400)
    diff = page.evaluate(
        "() => { const g = window.__viewer.gcPlot, m = window.__viewer.main;"
        " return [m.scales.x.min - g.scales.x.min, m.scales.x.max - g.scales.x.max, g.bbox.left - m.bbox.left, g.bbox.width - m.bbox.width]; }"
    )
    assert all(abs(v) <= 1 for v in diff)
    assert not page.errors


def test_page_without_optional_parts(open_page):
    page = open_page(with_gc=False, with_annotation=False, with_thresholds=False)
    assert page.locator("#gcplot canvas").count() == 0 and page.locator("#annplot canvas").count() == 0
    assert page.locator("#gcline").is_hidden()
    box = page.locator("#plot .u-over").bounding_box()
    page.mouse.move(box["x"] + 400, box["y"] + 100)
    page.mouse.wheel(0, -300)
    page.wait_for_timeout(400)
    assert view(page)[1] - view(page)[0] < LENGTH
    assert not page.errors
