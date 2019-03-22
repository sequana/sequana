# -*- coding: utf-8 -*-
#
#  This file is part of Sequana software
#
#  Copyright (c) 2019 - Sequana Development Team
#
#  File author(s):
#      Thomas Cokelaer <thomas.cokelaer@pasteur.fr>
#
#  Distributed under the terms of the 3-clause BSD license.
#  The full license is in the LICENSE file, distributed with this software.
#
#  website: https://github.com/sequana/sequana
#  documentation: http://sequana.readthedocs.io
#
##############################################################################
import json
from sequana.lazy import pylab
from sequana.lazy import pandas as pd

from sequana import logger
logger.name = __name__


__all__ = ["StatsFile"]


class StatsFile(object):
    """Reads a bcl2fastq Stats.json file


    """
    def __init__(self, filename="Stats.json"):
        with open(filename) as fh:
            self.data = json.load(fh)

    def to_summary_reads(self, filename):
        with open(filename, "w") as fh:
            fh.write("%4s %32s %16s\n" % ("Lane", "Sample", 'NumberReads'))
            for i, lane in enumerate(self.data["ConversionResults"]):
                total = 0
                for this in lane['DemuxResults']:
                    fh.write("%4s %32s %16s\n" %(i+1, this['SampleId'], this['NumberReads']))
                    total += this['NumberReads']
                fh.write("%4s %32s %16s\n" % (i+1, "Total", total))
                try:
                    fh.write("%4s %32s %16s\n" % (i+1, "Undetermined",
                        lane['Undetermined']['NumberReads']))
                except:
                    pass

    def get_data_reads(self):

        lanes = []
        names = []
        reads = []

        for i, lane in enumerate(self.data["ConversionResults"]):
            total = 0
            for this in lane['DemuxResults']:
                lanes.append(i+1)
                names.append(this['SampleId'])
                reads.append(this['NumberReads'])
                total += this['NumberReads']
            try:
                lanes.append(i+1)
                names.append("Undetermined")
                reads.append(lane['Undetermined']['NumberReads'])
            except:
                pass

        df = pd.DataFrame({"lane": lanes, "name": names, "count": reads})
        return df

    def barplot(self, filename="lane{}_status.png", lanes=None):
        df = self.get_data_reads()
        if lanes is None:
            lanes = df.lane.unique()

        for lane in lanes:
            pylab.clf()
            query = "lane==@lane and name!='Undetermined'"
            counts = df.query(query)['count']
            total = counts.sum()
            L = len(counts)

            query = "lane==@lane and name=='Undetermined'"
            under = df.query(query)['count'].sum()
            pylab.bar(range(L), counts, color="b", label="reads")

            if 100 * under / total <20:
                color = "green"
            elif 100 * under / total <50:
                color = "orange"
            else:
                color = "red"

            pylab.bar(range(L,L+1), under, color=color, label="undetermined")
            pylab.xticks([])
            pylab.ylabel("Number of reads")
            pylab.legend()
            pylab.title("Lane {}".format(lane))
            pylab.savefig(filename.format(lane), dpi=200)

    def barplot_summary(self, filename=None):
        df = self.get_data_reads()
        under = df.query("name=='Undetermined'")
        total = df.query("name!='Undetermined'")
        total = total.groupby("lane").sum().reset_index()
        total["name"] = "Total"
        df = pd.concat([under, total])
        df = df.pivot(index="lane", columns="name", values="count")
        df = df[["Total", "Undetermined"]]
        df.plot.barh(stacked=True)
        pylab.legend()
        pylab.xlabel("Number of reads")

        if filename:
            pylab.savefig(filename, dpi=200)
        return df




















