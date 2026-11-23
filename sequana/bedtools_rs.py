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
"""Rust-accelerated variants of the coverage classes.

The reference implementations (:class:`~sequana.bedtools.ChromosomeCov` and
:class:`~sequana.bedtools.FilteredGenomeCov`) are left untouched so that both
can be run side by side. Requires the ``sequana_cov_rs`` extension
(``rust/sequana_cov_rs``, install with ``maturin develop --release``).
"""
import copy
import re

import colorlog
import numpy as np
import pandas as pd

from sequana.bedtools import ChromosomeCov, FilteredGenomeCov, SequanaCoverage
from sequana.mixture_rs import EMRS

try:
    import sequana_cov_rs
except ImportError:  # pragma: no cover
    sequana_cov_rs = None


logger = colorlog.getLogger(__name__)


def rust_available():
    return sequana_cov_rs is not None


class FilteredGenomeCovRS(FilteredGenomeCov):
    """:class:`FilteredGenomeCov` whose ROI merging is done in Rust."""

    def _merge_region(self, zscore_label="zscore"):
        if sequana_cov_rs is None:  # pragma: no cover
            raise ImportError("sequana_cov_rs is not installed (see rust/sequana_cov_rs)")
        df = self.rawdf
        if len(df) == 0:
            return []
        out = sequana_cov_rs.merge_rois(
            df["pos"].to_numpy(dtype=np.int64),
            df["cov"].to_numpy(dtype=np.float64),
            df["rm"].to_numpy(dtype=np.float64),
            df[zscore_label].to_numpy(dtype=np.float64),
            float(self.thresholds.low),
            float(self.thresholds.high),
            int(self.step),
        )
        chrom = df["chr"].iloc[0]
        columns = {k: v.tolist() for k, v in out.items()}
        return [
            {
                "chr": chrom,
                "start": columns["start"][i],
                "end": columns["end"][i],
                "size": columns["size"][i],
                "mean_cov": columns["mean_cov"][i],
                "mean_rm": columns["mean_rm"][i],
                "mean_zscore": columns["mean_zscore"][i],
                "log2_ratio": columns["log2_ratio"][i],
                "max_zscore": columns["max_zscore"][i],
                "max_cov": columns["max_cov"][i],
            }
            for i in range(len(columns["start"]))
        ]


class ChromosomeCovRS(ChromosomeCov):
    """:class:`ChromosomeCov` using :class:`FilteredGenomeCovRS`."""

    _filtered_class = FilteredGenomeCovRS
    _em_class = EMRS

    def get_centralness(self, threshold=3):
        """Same as :meth:`ChromosomeCov.get_centralness`, computed in Rust.

        Centralness is 1 - (sum of ROI sizes) / length; the sign of a ROI does
        not matter so no dataframe of ROIs is built. Falls back to the Python
        implementation when annotations are present (a ROI may then be
        duplicated once per overlapping feature).
        """
        if sequana_cov_rs is None or self.bed.feature_dict or "zscore" not in self.df.columns:
            return super().get_centralness(threshold)

        df = self.df
        z = df["zscore"].to_numpy(dtype=np.float64)
        # the double thresholds follow low/high, as in the Python implementation
        th = copy.copy(self.thresholds)
        th.low, th.high = -threshold, threshold
        keep = (z > th.high2) | (z < th.low2)
        total = sequana_cov_rs.total_roi_size(
            np.ascontiguousarray(df["pos"].to_numpy(dtype=np.int64)[keep]),
            np.ascontiguousarray(z[keep]),
            float(-threshold),
            float(threshold),
            int(self.binning),
        )
        return 1 - total / float(len(self))

    _G_FORMAT = re.compile(r"^%\.(\d+)g$")

    def to_csv(self, filename=None, start=None, stop=None, **kwargs):
        """Same as :meth:`ChromosomeCov.to_csv`; the CSV string is built in Rust
        when it is requested in memory with ``index=False`` and a ``%.Ng``
        float format on numeric columns (used for the CanvasJS plots)."""
        fmt = self._G_FORMAT.match(str(kwargs.get("float_format", "")))
        if sequana_cov_rs is None or filename is not None or not fmt or kwargs.get("index", True) is not False:
            return super().to_csv(filename, start, stop, **kwargs)
        if set(kwargs) - {"columns", "index", "float_format"}:
            return super().to_csv(filename, start, stop, **kwargs)

        df = self.df.loc[start:stop]
        columns = list(kwargs.get("columns") or df.columns)
        arrays = []
        for name in columns:
            values = df[name]
            if values.dtype.kind in "iu":
                arrays.append(np.ascontiguousarray(values.to_numpy(dtype=np.int64)))
            elif values.dtype == np.float64:
                arrays.append(np.ascontiguousarray(values.to_numpy()))
            else:
                return super().to_csv(filename, start, stop, **kwargs)
        return sequana_cov_rs.to_csv([str(c) for c in columns], arrays, int(fmt.group(1)))

    def init(self):
        """Same as :meth:`ChromosomeCov.init`, but chunks come from the Rust BED reader."""
        offsets = getattr(self.bed, "_byte_offsets", None)
        if sequana_cov_rs is None or offsets is None or self.chrom_name not in offsets:
            return super().init()
        super().init()  # sets mode and metrics; the pandas iterator is lazy, replace it
        self._pandas_iterator = self.iterator
        self.iterator = self._rust_chunks(offsets[self.chrom_name], self.bed.positions[self.chrom_name]["N"])

    def _rust_chunks(self, offset, nrows):
        """Yield dataframes of at most *chunksize* rows; fall back to pandas if
        the file is not one the fast reader handles."""
        remaining = nrows
        first = True
        while remaining > 0:
            try:
                out = sequana_cov_rs.read_bed(self.bed.input_filename, offset, min(self.chunksize, remaining))
            except ValueError as err:
                if not first:
                    raise
                logger.warning(f"fast BED reader not used ({err}); falling back to pandas")
                yield from self._pandas_iterator
                return
            first = False
            n = len(out["pos"])
            if n == 0:
                return
            pos = out["pos"]
            columns = {"pos": pos, "cov": out["cov"].astype(np.int64) if out["cov_is_int"] else out["cov"]}
            if out["extra"] is not None:
                columns["mapq0"] = out["extra"].astype(np.int64) if out["extra_is_int"] else out["extra"]
            yield pd.DataFrame(columns, index=pd.Index(pos, name="pos"))
            offset = out["next_offset"]
            remaining -= n

    def _set_chunk(self, chunk):
        if 0 in chunk.columns:  # raw pandas chunk (fallback)
            return super()._set_chunk(chunk)
        # the chromosome column is not needed: the chunk is already indexed by position
        self._finish_chunk(chunk)


class SequanaCoverageRS(SequanaCoverage):
    """:class:`SequanaCoverage` scanning the BED file with the Rust reader."""

    _chromosome_class = ChromosomeCovRS

    def _scan_bed(self):
        if sequana_cov_rs is None:  # pragma: no cover
            return super()._scan_bed()
        logger.info("Scanning input file with the fast reader")
        try:
            names, first, last, offsets, total = sequana_cov_rs.scan_bed(self.input_filename)
        except Exception as err:  # pragma: no cover
            logger.warning(f"fast BED scan failed ({err}); falling back to pandas")
            return super()._scan_bed()
        self.chrom_names = list(names)
        self.total_length = total
        self.positions = {
            name: {"start": int(a), "end": int(b), "N": int(b) - int(a) + 1} for name, a, b in zip(names, first, last)
        }
        self._byte_offsets = dict(zip(names, offsets))
