#
#  This file is part of Sequana software
#
#  Copyright (c) 2016-2022 - Sequana Development Team
#
#  Distributed under the terms of the 3-clause BSD license.
#  The full license is in the LICENSE file, distributed with this software.
#
#  website: https://github.com/sequana/sequana
#  documentation: http://sequana.readthedocs.io
#
##############################################################################
import colorlog

from sequana.lazy import pandas as pd

logger = colorlog.getLogger(__name__)


__all__ = ["CheckM", "MultiCheckM"]


class CheckM:
    def __init__(self, filename):
        self.filename = filename
        # The format is not CSV or TSV, but a complex structure with different spaces...
        # The header can be just written here.
        header = [
            "sample",
            "marker_lineage",
            "#genomes",
            "#markers",
            "#marker_sets",
            "0",
            "1",
            "2",
            "3",
            "4",
            "5+",
            "Completeness",
            "Contamination",
            "Strain heterogeneity",
        ]

        with open(self.filename, "r") as fin:
            data = fin.read()
            values = data.split("\n")[3]

        # convert the string to number when possible
        new_values = []
        for val in values.split("  "):
            if val:
                try:
                    val = float(val)
                except ValueError:
                    try:
                        val = int(val)
                    except ValueError:
                        pass
                new_values.append(val)

        self.df = pd.Series(new_values)
        self.df.index = header


class MultiCheckM:
    def __init__(self, filenames):

        dfs = []

        for filename in filenames:
            try:
                dfs.append(CheckM(filename).df)
            except Exception:
                logger.warning(f"Skipped {filename}")

        self.df = pd.concat(dfs, axis=1)
