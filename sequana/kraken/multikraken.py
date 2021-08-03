#  This file is part of Sequana software
#
#  Copyright (c) 2016-2021 - Sequana Development Team
#
#  Distributed under the terms of the 3-clause BSD license.
#  The full license is in the LICENSE file, distributed with this software.
#
#  website: https://github.com/sequana/sequana
#  documentation: http://sequana.readthedocs.io
#
##############################################################################
import json

from sequana.lazy import pandas as pd
from sequana.lazy import pylab

import colorlog

logger = colorlog.getLogger(__name__)


__all__ = ["MultiKrakenResults", "MultiKrakenResults2"]


class MultiKrakenResults:
    """Select several kraken output and creates summary plots

    ::

        import glob
        mkr = MultiKrakenResults(glob.glob("*/*/kraken.csv"))
        mkr.plot_stacked_hist()

    """

    def __init__(self, filenames, sample_names=None):
        self.filenames = filenames
        if sample_names is None:
            self.sample_names = list(range(1, len(filenames) + 1))
        else:
            self.sample_names = sample_names

        self.ranks = [
            "kingdom",
            "phylum",
            "class",
            "order",
            "family",
            "genus",
            "species",
            "name",
        ]

    def get_df(self, limit=5):
        data = {}
        for sample, filename in zip(self.sample_names, self.filenames):
            df = pd.read_csv(filename)
            count = df["count"].sum()
            if "kingdom" not in df.columns:
                for name in self.ranks:
                    df[name] = "Unclassified"

            df = df.groupby("kingdom")["percentage"].sum()
            # if a taxon is obsolete, the kingdom is empty.
            # We will set the kingdom as Unclassified and raise a warning
            # if the count is > 5%
            if " " in df.index:
                percent = df.loc[" "]
                if percent > limit:
                    logger.warning(
                        "Found {}% of taxons with no kingdom".format(percent)
                    )
                if "Unclassified" in df.index:
                    df.loc["Unclassified"] += df.loc[" "]
                    df.drop(" ", inplace=True)
                else:
                    df.loc["Unclassified"] = df.loc[" "]
                    df.drop(" ", inplace=True)
            df["Count"] = count
            data[sample] = df

        df = pd.DataFrame(data)
        df = df.fillna(0)
        df = df.sort_index(ascending=False)
        df = df.sort_index(ascending=True, axis=1)

        return df

    def plot_stacked_hist(
        self,
        output_filename=None,
        dpi=200,
        kind="barh",
        fontsize=10,
        edgecolor="k",
        lw=1,
        width=1,
        ytick_fontsize=10,
        max_labels=50,
    ):

        """Summary plot of reads classified."""
        df = self.get_df()

        # remove count from the summary graph
        df.drop("Count", inplace=True)

        fig, ax = pylab.subplots(figsize=(9.5, 7))

        labels = []
        for kingdom in sorted(df.index):
            if kingdom == "Eukaryota":
                color = "purple"
            elif kingdom == "Unclassified":
                color = "grey"
            elif kingdom == "Bacteria":
                color = "red"
            elif kingdom == "Viruses":
                color = "darkgreen"
            elif kingdom == "Archaea":
                color = "yellow"
            else:
                color = "darkblue"
            if kind == "barh":
                pylab.barh(
                    range(0, len(df.columns)),
                    df.loc[kingdom],
                    height=width,
                    left=df.loc[labels].sum().values,
                    edgecolor=edgecolor,
                    lw=lw,
                    color=color,
                    alpha=0.8,
                )
            else:
                pylab.bar(
                    range(0, len(df.columns)),
                    df.loc[kingdom],
                    width=width,
                    bottom=df.loc[labels].sum().values,
                    edgecolor=edgecolor,
                    lw=lw,
                    color=color,
                    alpha=0.8,
                )
            labels.append(kingdom)

        if kind == "barh":
            pylab.xlabel("Percentage (%)", fontsize=fontsize)
            pylab.ylabel("Sample index/name", fontsize=fontsize)
            if len(self.sample_names) < max_labels:
                pylab.yticks(
                    range(len(self.sample_names)),
                    self.sample_names[::-1],
                    fontsize=ytick_fontsize,
                )
            else:
                pylab.yticks([1], [""])
            pylab.xlim([0, 100])
            pylab.ylim([-0.5, len(df.columns) - 0.5])
        else:
            pylab.ylabel("Percentage (%)", fontsize=fontsize)
            pylab.xlabel("Sample index/name", fontsize=fontsize)
            pylab.ylim([0, 100])
            if len(self.sample_names) < max_labels:
                pylab.xticks(
                    range(len(self.sample_names)),
                    self.sample_names[::-1],
                    fontsize=ytick_fontsize,
                )
            else:
                pylab.xticks([1], [""])
            pylab.xlim([-0.5, len(df.columns) - 0.5])

        ax.legend(labels, title="kingdom", bbox_to_anchor=(1, 1))

        try:  # pragma: no cover
            pylab.tight_layout()
        except:  # pragma: no cover
            pass

        if output_filename:
            pylab.savefig(output_filename, dpi=dpi)


class MultiKrakenResults2:
    """Select several kraken output and creates summary plots

    ::

        import glob
        mkr = MultiKrakenResults2(glob.glob("*/*/summary.json"))
        mkr.plot_stacked_hist()

    """

    def __init__(self, filenames, sample_names=None):

        self.filenames = filenames
        if sample_names is None:
            self.sample_names = list(range(1, len(filenames) + 1))
        else:
            self.sample_names = sample_names

    def get_df(self, limit=5, sorting_method="sample_name"):
        data = {}
        for sample, filename in zip(self.sample_names, self.filenames):
            summary = json.loads(open(filename, "r").read())
            total = summary["total"]
            if "unclassified" not in summary:
                summary["unclassified"] = 0
            data[sample] = {
                "unclassified": round(summary["unclassified"] / total * 100, 2),
                "nreads": summary["total"],
            }
            for db in summary["databases"]:
                try:
                    data[sample][db] = round(summary[db]["C"] / total * 100, 2)
                except:
                    data[sample][db] = 0
        df = pd.DataFrame(data)
        df = df.fillna(0)
        df = df.loc[["unclassified"] + [x for x in df.index if x != "unclassified"]]
        # df = df.sort_index(ascending=False)
        df = df.sort_index(ascending=False, axis=1)
        return df

    def plot_stacked_hist(
        self,
        output_filename=None,
        dpi=200,
        kind="barh",
        fontsize=10,
        edgecolor="k",
        lw=2,
        width=1,
        ytick_fontsize=10,
        max_labels=50,
        logx=False,
        alpha=0.8,
        colors=None,
        cmap="hot_r",
        sorting_method="sample_name",
    ):
        """Summary plot of reads classified.


        :param sorting_method: only by sample name for now
        """
        df = self.get_df(sorting_method=sorting_method)
        df = df.T
        del df["nreads"]

        fig, ax = pylab.subplots(figsize=(9.5, 7))

        # we will store grey/unclassified first and then other DB with a max of
        # 10 DBs
        L = len(df.columns) - 1
        from matplotlib import cm

        if colors is None:
            colors = [cm.get_cmap(cmap)(x) for x in pylab.linspace(0.2, 1, L)]
        colors = ["grey"] + colors
        df.plot(
            kind="barh",
            stacked=True,
            width=width,
            edgecolor=edgecolor,
            color=colors,
            lw=lw,
            alpha=alpha,
            ax=ax,
            legend=False,
        )
        if logx is True:
            pylab.semilogx()

        pylab.xlabel("Percentage (%)", fontsize=fontsize)
        pylab.ylabel("Sample index/name", fontsize=fontsize)
        if len(self.sample_names) < max_labels:
            pylab.yticks(
                range(len(self.sample_names)), df.index, fontsize=ytick_fontsize
            )
        else:
            pylab.yticks([1], [""])
        pylab.xlim([0, 100])
        pylab.ylim([0 - 0.5, len(df.index) - 0.5])  # +1 for the legends
        ax.legend(df.columns, title="DBs ", bbox_to_anchor=(1, 1))

        try:  # pragma: no cover
            pylab.tight_layout()
        except:  # pragma: no cover
            pass

        if output_filename:
            pylab.savefig(output_filename, dpi=dpi)
        return df
