#
#  This file is part of Sequana software
#
#  Copyright (c) 2016-2021 - Sequana Development Team
#
#  File author(s): Sequana team
#
#  Distributed under the terms of the 3-clause BSD license.
#  The full license is in the LICENSE file, distributed with this software.
#
#  website: https://github.com/sequana/sequana
#  documentation: http://sequana.readthedocs.io
#
##############################################################################
from sequana.utils.pandas import PandasReader

from sequana.lazy import pylab

__all__ = ["FRiP"]


class FRiP:
    """Reader for output of FRiP computation in chipseq pipeline

    FRiP is the fraction of reads found in Peaks.

    expected format::

        bamfile,count,in_peaks, FRiP,
        1_S1_L001.sorted.bam,3548090,53673,0.015127293839784221
        2_S2_L001.sorted.bam,3868608,58292,0.01506795209026089
        3_S3_L001.sorted.bam,4092990,50219,0.01226951446253228

    """

    def __init__(self, filename):
        self.df = PandasReader(filename, sep=",").df

    def plot(self):
        """scatter plot of FRiP verus Reads in peaks """

        pylab.clf()
        MX = self.df.FRiP.max()
        MY = self.df["in_peaks"].max()
        pylab.plot([0, MX], [0, MY], ls="--", color="b", alpha=0.5)
        for bamfile in self.df["bamfile"]:
            label = bamfile.split("/")[0]
            self.df.query("bamfile==@bamfile").plot(
                x="FRiP", y="in_peaks", marker="o", alpha=0.5, lw=0, label=label, ax=pylab.gca()
            )
        pylab.ylabel("Reads in peaks")
        pylab.xlabel("FRiP")
        pylab.xlim(0, pylab.xlim()[1])
        pylab.ylim(0, pylab.ylim()[1])
        try:
            pylab.tight_layout()
        except Exception: #pragma: no cover
            pass
        pylab.grid()

