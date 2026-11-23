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
"""Interactive coverage page for one chromosome analysed by ``sequana_coverage``.

This is the only module of :mod:`sequana.viewers` that knows about the coverage
analysis. It reads a :class:`sequana.bedtools.ChromosomeCov` (duck-typed: ``df``,
``thresholds``, ``best_gaussian``, ``DOC``, ``BOC``) and hands plain arrays to
:func:`sequana.viewers.genome_track.build_page`.
"""
import os

import colorlog
import numpy as np

from sequana.viewers.genome_track import build_page, load_annotation

logger = colorlog.getLogger(__name__)

__all__ = ["coverage_page", "write_coverage_viewer"]


def coverage_page(chrom, results, window, annotation_file=None, command=""):
    """HTML page of one chromosome, or None when the chromosome cannot be shown.

    :param chrom: :class:`~sequana.bedtools.ChromosomeCov` after ``run()``
    :param results: what ``chrom.run()`` returned (gives the regions of interest)
    :param int window: running median window
    :param annotation_file: optional GFF3 or GenBank file (gene track and gene names)
    """
    name = chrom.chrom_name
    if getattr(chrom, "_mode", "memory") != "memory":
        logger.warning(f"Viewer skipped for {name}: the chromosome is read by chunks (above --chunk-size)")
        return None
    df = chrom.df
    pos = df["pos"].to_numpy()
    if len(pos) == 0 or not (np.diff(pos) == 1).all():
        logger.warning(f"Viewer skipped for {name}: positions are not contiguous (is --binning used?)")
        return None

    # NaN at both ends (median and GC windows): extend the nearest value
    baseline = df["rm"].ffill().bfill().to_numpy(dtype=np.float64)
    gc = df["gc"].ffill().bfill().to_numpy(dtype=np.float64) if "gc" in df.columns else None

    table = results.get_rois().df
    # with an annotation, sequana repeats a region once per overlapping feature
    table = table.drop_duplicates(subset=["start", "end"])
    regions = [
        {
            "start": int(r.start),
            "end": int(r.end),
            "size": int(r["size"]),
            "mean_signal": float(r.mean_cov),
            "mean_baseline": float(r.mean_rm),
            "log2_ratio": float(r.log2_ratio) if r.log2_ratio is not None else None,
            "score": float(r.max_zscore),
            "type": "low" if r.max_zscore < 0 else "high",
        }
        for _, r in table.iterrows()
    ]

    features, lanes = None, None
    if annotation_file:
        features, lanes = load_annotation(annotation_file, name)
        if not features:
            logger.warning(f"No feature found for {name} in {annotation_file}: no gene track")
            features = None

    gaussian = chrom.best_gaussian
    thresholds = {
        "low": float(chrom.thresholds.low),
        "high": float(chrom.thresholds.high),
        "mu": float(gaussian["mu"]),
        "sigma": float(gaussian["sigma"]),
    }
    stats = [
        ("length", f"{len(pos):,} bp"),
        ("mean depth", f"{chrom.DOC:.1f}"),
        ("breadth", f"{chrom.BOC:.1f}%"),
        ("median window", f"{window:,}"),
        ("thresholds", f"[{thresholds['low']:g}, {thresholds['high']:g}]"),
    ]
    return build_page(
        name,
        df["cov"].to_numpy(dtype=np.float64),
        baseline,
        regions,
        pos0=int(pos[0]),
        gc=gc,
        features=features,
        lanes=lanes,
        thresholds=thresholds,
        signal_mean=float(chrom.DOC),
        stats=stats,
        labels={"title": "Coverage analysis of chromosome {chrom}"},
        command=command,
    )


def write_coverage_viewer(chrom, results, directory, window, annotation_file=None, command=""):
    """Write ``<directory>/<chromosome>.viewer.html``.

    :return: the file name (relative to *directory*), or None if the page was not created
    """
    page = coverage_page(chrom, results, window, annotation_file=annotation_file, command=command)
    if page is None:
        return None
    filename = f"{chrom.chrom_name}.viewer.html"
    with open(os.path.join(directory, filename), "w") as fout:
        fout.write(page)
    logger.info(f"Interactive viewer written in {os.path.join(directory, filename)}")
    return filename
