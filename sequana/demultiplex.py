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
    """Reads a bcl2fastq Stats.json file and produces useful plots for QCs
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
                lanes.append(lane['LaneNumber'])
                names.append(this['SampleId'])
                reads.append(this['NumberReads'])
                total += this['NumberReads']

            # if only undetermined (no sample sheet), no data should be found
            # meaning that there is 0 determined and no names associated
            # so we call it determined and store 0
            if total == 0:
                assert total ==0
                names.append("Determined")
                reads.append(0)
                lanes.append(lane["LaneNumber"])

            if "Undetermined" in lane:
                names.append("Undetermined")
                lanes.append(i+1)
                reads.append(lane['Undetermined']['NumberReads'])
            else:
                names.append("Undetermined")
                lanes.append(i+1)
                reads.append(0)
        #print(lanes, names, reads)

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
            if total >0:
                pylab.bar(range(L), counts, color="b", label="reads")

            if total == 0:
                color = "red"
            else:
                if 100 * under / total <20:
                    color = "green"
                elif 100 * under / total <50:
                    color = "orange"
                else:
                    color = "red"

            pylab.bar(range(L,L+1), under, color=color, label="undetermined")
            pylab.xticks([])
            pylab.ylabel("Number of reads")
            try:pylab.legend(loc="best")
            except:pass
            pylab.title("Lane {}".format(lane))
            pylab.savefig(filename.format(lane), dpi=200)

    def barplot_summary(self, filename=None, color=["green", "red"], 
            alpha=0.8):

        df = self.get_data_reads()
        under = df.query("name=='Undetermined'")

        total = df.query("name!='Undetermined'")
        total = total.groupby("lane").sum().reset_index()
        total["name"] = "Determined"

        df = pd.concat([under, total], sort=True)

        df = df.pivot(index="lane", columns="name", values="count")
        df = df[["Determined", "Undetermined"]]
        df.plot.barh(stacked=True, color=color, alpha=alpha)
        pylab.legend()
        pylab.xlabel("Number of reads")

        if filename:
            pylab.savefig(filename, dpi=200)
        return df

    def plot_unknown_barcodes(self, N=20):
        ub = self.data['UnknownBarcodes']
        df = pd.DataFrame({x['Lane']:x['Barcodes'] for x in ub})
        if "unknown" in df.index and len(df) == 1:
            df.loc['known'] = [0 for i in df.columns]

        # if data is made of undetermined only, the dataframe is just made of
        # N lanes with one entry : unknown
        S = df.sum(axis=1).sort_values(ascending=False).index[0:N]
        data = df.loc[S][::-1]
        #print(data)

        data.columns = ["Lane {}".format(x) for x in data.columns]
        from matplotlib import rcParams
        rcParams['axes.axisbelow'] = True
        pylab.figure(figsize=(10,8))
        ax = pylab.gca()
        data.plot(kind="barh", width=1, ec="k", ax=ax)
        rcParams['axes.axisbelow'] = False
        pylab.xlabel("Number of reads", fontsize=12)
        pylab.ylabel("Frequency (%)", fontsize=12)
        pylab.grid(True)
        pylab.legend(["Lane {}".format(x) for x in range(1, len(df.columns)+1)])
        try:
            pylab.tight_layout()
        except Exception as err:
            print(err)
        return data
