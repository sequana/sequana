import json
import re

import numpy as np
import pytest

from sequana.viewers.genome_track import (
    annotate_regions,
    build_page,
    build_pyramid,
    decode_array,
    encode_array,
    load_annotation,
    write_page,
)

GFF = """##gff-version 3
chrA\tsrc\tregion\t1\t5000\t.\t.\t.\tID=chrA
chrA\tsrc\tgene\t100\t900\t.\t+\t.\tID=gene1;Name=alpha;product=first%2C protein
chrA\tsrc\tmRNA\t100\t900\t.\t+\t.\tID=gene1.mRNA;Parent=gene1
chrA\tsrc\tCDS\t100\t900\t.\t+\t.\tID=gene1.CDS;Parent=gene1.mRNA
chrA\tsrc\tgene\t500\t1500\t.\t-\t.\tID=gene2;gene_id=beta
chrA\tsrc\ttRNA\t2000\t2075\t.\t+\t.\tID=trna1;product=tRNA-Ala
chrB\tsrc\tgene\t10\t90\t.\t+\t.\tID=other
"""


@pytest.fixture
def gff(tmp_path):
    path = tmp_path / "test.gff3"
    path.write_text(GFF)
    return path


def regions():
    return [
        {"start": 120, "end": 130, "mean_signal": 5.0, "mean_baseline": 30.0, "score": -8.0, "type": "low"},
        {
            "start": 1600,
            "end": 1700,
            "mean_signal": 90.0,
            "mean_baseline": 30.0,
            "score": 9.0,
            "type": "high",
            "log2_ratio": 1.58,
        },
        {"start": 3000, "end": 3001, "mean_signal": 5.0, "mean_baseline": 30.0, "score": -5.0, "type": "low"},
    ]


def test_encode_decode_roundtrip():
    rng = np.random.default_rng(0)
    values = np.concatenate([rng.poisson(30, 1000).astype(float), [0.0, 1e6, 0.25]])
    decoded = decode_array(encode_array(values), len(values))
    np.testing.assert_array_equal(decoded, values.astype(np.float32))


def test_pyramid_levels_summarise_the_signal():
    rng = np.random.default_rng(1)
    signal = rng.poisson(30, 100003).astype(float)
    baseline = np.full(len(signal), 30.0)
    gc = np.linspace(0.3, 0.7, len(signal))
    levels = build_pyramid(signal, baseline, gc)

    assert levels[0]["bin"] == 1 and set(levels[0]["arrays"]) == {"cov", "rm"}
    assert [lv["bin"] for lv in levels] == [4**i for i in range(len(levels))]
    assert levels[-1]["n"] <= 1500 * 4  # coarsest level is small
    for level in levels[1:]:
        n, size = level["n"], level["bin"]
        assert n == int(np.ceil(len(signal) / size))
        arrays = {k: decode_array(v, n) for k, v in level["arrays"].items()}
        first = signal[:size]  # first bin
        assert arrays["min"][0] == first.min() and arrays["max"][0] == first.max()
        assert arrays["mean"][0] == pytest.approx(first.mean(), rel=1e-6)
        assert ("gc" in arrays) == (size >= 16)
        if "gc" in arrays:  # quantised to 1/256
            assert arrays["gc"][0] == pytest.approx(gc[:size].mean(), abs=1 / 256)
    # peaks survive the reduction (the point of min/max bins)
    signal[5000] = 12345
    levels = build_pyramid(signal, baseline)
    top = decode_array(levels[-1]["arrays"]["max"], levels[-1]["n"])
    assert top.max() == 12345


def test_load_annotation_gff(gff):
    features, lanes = load_annotation(gff, "chrA")
    # region, mRNA and CDS are skipped; genes and non-coding RNAs are kept
    assert [f["id"] for f in features] == ["gene1", "gene2", "trna1"]
    assert features[0]["name"] == "alpha" and features[0]["product"] == "first, protein"  # %2C decoded
    assert features[1]["name"] == "beta"  # falls back to gene_id
    assert features[2]["type"] == "tRNA"
    # gene1 and gene2 overlap: two lanes
    assert lanes == 2 and features[0]["lane"] != features[1]["lane"] and features[2]["lane"] == 0
    assert load_annotation(gff, "chrB")[0][0]["id"] == "other"
    assert load_annotation(gff, "missing") == ([], 1)


def test_load_annotation_falls_back_to_cds(tmp_path):
    path = tmp_path / "cds_only.gff3"
    path.write_text("chrA\tsrc\tCDS\t10\t99\t.\t+\t.\tID=c1;Name=x\nchrA\tsrc\tCDS\t200\t299\t.\t+\t.\tID=c2\n")
    features, _ = load_annotation(path, "chrA")
    assert [f["name"] for f in features] == ["x", "c2"]


def test_annotate_regions_uses_exclusive_end(gff):
    features, _ = load_annotation(gff, "chrA")
    hits = annotate_regions(regions(), features)
    assert hits[0] == [0]  # 120-130 is inside gene1 only
    assert hits[1] == []  # 1600-1700 does not overlap gene2 (500-1500)
    assert hits[2] == []  # intergenic
    # the region end is exclusive: [50, 100) stops one base before gene1 (starts at 100)
    edge = [
        {"start": 50, "end": 100},
        {"start": 50, "end": 101},
        {"start": 900, "end": 950},
        {"start": 901, "end": 950},
    ]
    assert annotate_regions(edge, features) == [[], [0], [0, 1], [1]]
    assert annotate_regions(regions(), []) == [[], [], []]


def _meta(page):
    text = re.search(r'<script id="meta" type="application/json">(.*?)</script>', page, re.S).group(1)
    return json.loads(text.replace("<\\/", "</"))


def test_build_page_content(gff):
    rng = np.random.default_rng(2)
    signal = rng.poisson(30, 5000).astype(float)
    features, lanes = load_annotation(gff, "chrA")
    page = build_page(
        "chrA",
        signal,
        np.full(5000, 30.0),
        regions(),
        gc=np.full(5000, 0.5),
        features=features,
        lanes=lanes,
        thresholds={"low": -4, "high": 4, "mu": 1.0, "sigma": 0.1},
        stats=[("length", "5,000 bp")],
        command="tool <-x> & more",
    )
    assert page.startswith("<!DOCTYPE html>")
    assert not re.search(r"__[A-Z_]+__", re.sub(r"<script id=\"meta\".*?</script>", "", page, flags=re.S)[:2000])
    meta = _meta(page)
    assert meta["name"] == "chrA" and meta["length"] == 5000 and meta["pos0"] == 1
    assert meta["has_gc"] and meta["thresholds"]["sigma"] == 0.1
    assert len(meta["ann"]["features"]) == 3 and meta["ann"]["lanes"] == 2
    assert [r[7] for r in meta["rois"]] == ["low", "high", "low"]
    assert meta["rois"][0][8] == [0] and meta["rois"][2][8] == []  # gene indices
    assert meta["rois"][0][5] is None and meta["rois"][1][5] == pytest.approx(1.58)  # log2 ratio
    assert "tool &lt;-x&gt; &amp; more" in page  # the command is escaped
    assert "3 regions" not in page and "(2 below, 1 above)" in page
    # without optional parts
    bare = _meta(build_page("chrA", signal, np.full(5000, 30.0), regions()))
    assert bare["thresholds"] is None and bare["ann"] is None and not bare["has_gc"]
    assert all(len(r) == 8 for r in bare["rois"])


def test_build_page_rejects_bad_input():
    signal = np.ones(100)
    with pytest.raises(ValueError):
        build_page("x", signal, np.ones(99), [])
    with pytest.raises(ValueError):
        build_page("x", signal, np.r_[np.nan, np.ones(99)], [])
    with pytest.raises(ValueError):
        build_page("x", signal, np.ones(100), [], gc=np.ones(50))


def test_script_end_tag_in_data_cannot_break_the_page(tmp_path):
    features = [
        {
            "start": 1,
            "end": 5,
            "strand": "+",
            "type": "gene",
            "name": "</script><b>",
            "product": "",
            "id": "x",
            "lane": 0,
        }
    ]
    path = write_page(tmp_path / "sub" / "page.html", "c", np.ones(50), np.ones(50), [], features=features, lanes=1)
    page = path.read_text()
    assert "</script><b>" not in page and _meta(page)["ann"]["features"][0][5] == "</script><b>"
