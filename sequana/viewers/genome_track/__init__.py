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
"""Generic multi-resolution viewer for a per-base signal (coverage, depth, ...).

Give it numpy arrays, a table of regions and optionally an annotation, get one
self-contained HTML page (no server, works from ``file://``)::

    from sequana.viewers.genome_track import build_page, load_annotation

    features, lanes = load_annotation("genes.gff3", "chr1")
    html = build_page("chr1", signal, baseline, regions, features=features, lanes=lanes)

The page uses uPlot (MIT, vendored in ``assets/uplot``).
"""
from .annotation import annotate_regions, load_annotation
from .builder import build_page, build_pyramid, decode_array, encode_array, write_page

__all__ = [
    "annotate_regions",
    "build_page",
    "build_pyramid",
    "decode_array",
    "encode_array",
    "load_annotation",
    "write_page",
]
