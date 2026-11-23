#
#  This file is part of Sequana software
#
#  Copyright (c) 2016-2026 - Sequana Development Team
#
#  Distributed under the terms of the 3-clause BSD license.
#  The full license is in the LICENSE file, distributed with this software.
#
#  website: https://github.com/sequana/sequana
#  documentation: http://sequana.readthedocs.io
#
##############################################################################
"""Build the self-contained HTML page of the genome-track viewer.

The per-base arrays are reduced into a pyramid of (min, max, mean) bins, byte-shuffled,
deflated and embedded in the page as base64; the JavaScript (uPlot) picks the level
that matches the zoom.
"""
import base64
import html
import json
import re
import zlib
from pathlib import Path

import numpy as np

__all__ = ["build_page", "build_pyramid", "decode_array", "encode_array", "write_page"]

ASSETS = Path(__file__).parent / "assets"
FACTOR = 4  # bin size ratio between two levels
MIN_BINS = 1500  # stop building coarser levels below this many bins
GC_MIN_BIN = 16  # the GC curve is smooth: only stored from this bin size on
MAX_LANES = 6

DEFAULT_LABELS = {
    "signal": "Coverage",
    "y_axis": "Per-base coverage",
    "baseline": "Running median",
    "gc": "GC content",
    "mean_signal": "Mean cov.",
    "mean_baseline": "Mean median",
    "score": "Max z-score",
    "low": "below average",
    "high": "above average",
    "low_short": "below",
    "high_short": "above",
    "title": "Coverage analysis of {chrom}",
}


def encode_array(array):
    """float array -> base64(raw deflate(byte-shuffled float32)).

    Shuffling groups the bytes of the same significance together, which compresses much
    better. :func:`decode_array` and the JavaScript do the reverse.
    """
    array = np.ascontiguousarray(array, dtype="<f4")
    shuffled = array.view(np.uint8).reshape(-1, 4).T.tobytes()
    compressor = zlib.compressobj(9, zlib.DEFLATED, -15)
    return base64.b64encode(compressor.compress(shuffled) + compressor.flush()).decode("ascii")


def decode_array(text, n):
    """Inverse of :func:`encode_array` (*n* is the number of values)."""
    shuffled = np.frombuffer(zlib.decompress(base64.b64decode(text), -15), dtype=np.uint8)
    return shuffled.reshape(4, n).T.copy().view("<f4").reshape(n)


def _quantise_gc(gc):
    """GC fraction rounded to 1/256 (exact in float32, so it compresses well)."""
    return np.round(gc * 256) / 256


def build_pyramid(signal, baseline, gc=None, factor=FACTOR, min_bins=MIN_BINS, gc_min_bin=GC_MIN_BIN):
    """Multi-resolution representation of the arrays.

    :return: list of levels. Level 0 has ``bin == 1`` and the arrays ``cov`` and ``rm``
        (and nothing else: min = max = mean = signal). The other levels hold the ``min``,
        ``max`` and ``mean`` of the signal, the mean ``rm`` and, from *gc_min_bin* on, the
        mean ``gc``. Arrays are encoded with :func:`encode_array`.
    """
    signal = np.asarray(signal, dtype=np.float64)
    baseline = np.asarray(baseline, dtype=np.float64)
    n = len(signal)
    levels = [{"bin": 1, "n": n, "arrays": {"cov": encode_array(signal), "rm": encode_array(baseline)}}]
    size = 1
    while n / size > min_bins:
        size *= factor
        edges = np.arange(0, n, size)
        counts = np.diff(np.append(edges, n))
        arrays = {
            "min": np.minimum.reduceat(signal, edges),
            "max": np.maximum.reduceat(signal, edges),
            "mean": np.add.reduceat(signal, edges, dtype=np.float64) / counts,
            "rm": np.add.reduceat(baseline, edges, dtype=np.float64) / counts,
        }
        if gc is not None and size >= gc_min_bin:
            arrays["gc"] = _quantise_gc(np.add.reduceat(gc, edges, dtype=np.float64) / counts)
        levels.append({"bin": size, "n": len(edges), "arrays": {k: encode_array(v) for k, v in arrays.items()}})
    if gc is not None and not any("gc" in level["arrays"] for level in levels):
        # short sequence: no coarse level was built, keep the GC at base resolution
        levels[0]["arrays"]["gc"] = encode_array(_quantise_gc(gc))
    return levels


def _region_rows(regions, gene_hits):
    """[start, end, size, mean_signal, mean_baseline, log2, score, type, genes]"""
    rows = []
    for i, region in enumerate(regions):
        log2 = region.get("log2_ratio")
        rows.append(
            [
                int(region["start"]),
                int(region["end"]),
                int(region.get("size", region["end"] - region["start"])),
                float(region["mean_signal"]),
                float(region["mean_baseline"]),
                None if log2 is None or not np.isfinite(log2) else float(log2),
                float(region["score"]),
                region["type"],
            ]
            + ([gene_hits[i]] if gene_hits is not None else [])
        )
    return rows


def build_page(
    name,
    signal,
    baseline,
    regions,
    *,
    pos0=1,
    gc=None,
    features=None,
    lanes=None,
    thresholds=None,
    signal_mean=None,
    stats=None,
    labels=None,
    command="",
):
    """Build the HTML page of one sequence.

    :param str name: sequence (chromosome) name
    :param signal: per-base signal (e.g. coverage), one value per position, no gaps
    :param baseline: per-base baseline (e.g. running median), same length, no NaN
    :param regions: list of dicts with ``start``, ``end`` (exclusive), ``mean_signal``,
        ``mean_baseline``, ``score``, ``type`` (``"low"`` or ``"high"``) and optionally
        ``size`` and ``log2_ratio``
    :param int pos0: position (1-based) of the first value
    :param gc: optional GC fraction (0-1) per position, no NaN: drawn on a right axis
    :param features: optional annotation, as returned by
        :func:`~sequana.viewers.genome_track.load_annotation`; gives a gene track and
        gene names in the region table
    :param lanes: number of lanes of *features* (from ``load_annotation``)
    :param thresholds: optional dict ``low``, ``high``, ``mu``, ``sigma``: the confidence
        band drawn is ``(threshold * sigma + mu) * baseline``
    :param signal_mean: mean level of the signal (used to clip the y axis); default: computed
    :param stats: list of ``(label, value)`` shown in the page header
    :param labels: overrides of :data:`DEFAULT_LABELS` (axis, legend and table texts)
    :param command: text of the "Command" section (hidden when empty)
    :return: the page as a string
    """
    from .annotation import annotate_regions

    signal = np.asarray(signal, dtype=np.float64)
    baseline = np.asarray(baseline, dtype=np.float64)
    if len(signal) != len(baseline) or (gc is not None and len(gc) != len(signal)):
        raise ValueError("signal, baseline and gc must have the same length")
    if np.isnan(baseline).any() or (gc is not None and np.isnan(gc).any()):
        raise ValueError("baseline and gc must not contain NaN (fill the edges first)")
    labels = dict(DEFAULT_LABELS, **(labels or {}))

    ann = None
    gene_hits = None
    if features:
        gene_hits = annotate_regions(regions, features)
        ann = {
            "lanes": min(lanes or 1, MAX_LANES),
            "features": [
                [
                    f["start"],
                    f["end"],
                    f["strand"],
                    min(f["lane"], MAX_LANES - 1),
                    f["type"],
                    f["name"],
                    f["product"],
                    f["id"],
                ]
                for f in features
            ],
        }
    rows = _region_rows(regions, gene_hits)
    meta = {
        "name": name,
        "pos0": int(pos0),
        "length": len(signal),
        "thresholds": thresholds,
        "signal_mean": float(signal_mean if signal_mean is not None else np.mean(signal)),
        "labels": labels,
        "has_gc": gc is not None,
        "ann": ann,
        "levels": build_pyramid(signal, baseline, None if gc is None else np.asarray(gc, dtype=np.float64)),
        "rois": rows,
    }
    n_low = sum(r[7] == "low" for r in rows)
    header = list(stats or [])
    if ann:
        header.append(("features", f"{len(ann['features']):,}"))
    header.append(
        ("regions", f"{len(rows)} ({n_low} {labels['low_short']}, {len(rows) - n_low} {labels['high_short']})")
    )

    page = (ASSETS / "template.html").read_text()
    substitutions = {
        "__TITLE__": html.escape(labels["title"].format(chrom=name)),
        "__STATS__": "".join(f"<span>{html.escape(str(k))} <b>{html.escape(str(v))}</b></span>" for k, v in header),
        "__COMMAND_SECTION__": f"<h2>Command</h2>\n  <pre>{html.escape(command)}</pre>" if command else "",
        "__ROI_LOW__": html.escape(labels["low"]),
        "__ROI_HIGH__": html.escape(labels["high"]),
        "__UPLOT_CSS__": (ASSETS / "uplot" / "uPlot.min.css").read_text(),
        "__VIEWER_CSS__": (ASSETS / "viewer.css").read_text(),
        "__UPLOT_JS__": (ASSETS / "uplot" / "uPlot.iife.min.js").read_text(),
        "__VIEWER_JS__": (ASSETS / "viewer.js").read_text(),
        # </ would end the script element
        "__META__": json.dumps(meta, separators=(",", ":")).replace("</", "<\\/"),
    }
    # one pass: inserted content (JS, JSON) is never scanned for placeholders
    return re.sub(r"__[A-Z_]+__", lambda m: substitutions.get(m.group(0), m.group(0)), page)


def write_page(path, *args, **kwargs):
    """:func:`build_page`, written to *path* (parent directories are created)."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(build_page(*args, **kwargs))
    return path
