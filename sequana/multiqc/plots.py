#
#  This file is part of Sequana software
#
#  Copyright (c) 2016-2021 - Sequana Dev Team (https://sequana.readthedocs.io)
#
#  Distributed under the terms of the 3-clause BSD license.
#  The full license is in the LICENSE file, distributed with this software.
#
#  Website:       https://github.com/sequana/sequana
#  Documentation: http://sequana.readthedocs.io
#  Contributors:  https://github.com/sequana/sequana/graphs/contributors
##############################################################################
import pandas as pd


class Reader:
    def __init__(self, filename):
        self.filename = filename
        self.reader()

    def reader(self):
        self.df = pd.read_csv(self.filename, sep="\t")


class Bowtie1Reader(Reader):
    def __init__(self, filename):
        super().__init__(filename)

    def plot_bar(self, html_code=False):
        import plotly.graph_objects as go

        fig = go.Figure()

        x = self.df["not_aligned_percentage"]
        compx = 100 - x

        fig.add_trace(go.Bar(y=self.df["Sample"], x=compx, name="Aligned", marker_color="#389f1f", orientation="h"))
        fig.add_trace(
            go.Bar(
                y=self.df.Sample,
                x=self.df["not_aligned_percentage"],
                name="Not Aligned",
                marker_color="#120946",
                orientation="h",
            )
        )

        # Here we modify the tickangle of the xaxis, resulting in rotated labels.
        fig.update_layout(barmode="stack", xaxis_tickangle=-45, height=400, title="Mapping on ribosomal/contaminant")

        if html_code:
            return fig
        else:  # pragma: no cover
            fig.show()


class Bowtie2(Reader):
    def __init__(self, filename):
        super().__init__(filename)

    def plot(self, html_code=False):
        import plotly.graph_objects as go

        fig = go.Figure()

        # get percentage instead of counts
        if "unpaired_aligned_multi" in self.df:
            # version <=0.17.2 we used the mqc_bowtie output
            # in version >0.17.2  we use multiqc_bowtie2.txt
            # self.df[["SE mapped uniquely", "SE multimapped", "SE not aligned"]].sum(axis=1).values
            # then printed in this order: SE mapped uniquely, SE multimapped, SE not aligned

            columns = [
                "unpaired_aligned_one",
                "unpaired_aligned_multi",
                "unpaired_aligned_none",
            ]
            names = ["Mapped", "Multi mapped", "Unmapped"]

            colors = ["#17478f", "#e2780d", "#9f1416"]

            S = self.df[columns].sum(axis=1).values
            df = self.df[columns].divide(S, axis=0) * 100
            df["Sample"] = self.df["Sample"]

            fig = go.Figure(
                data=[
                    go.Bar(
                        name=name,
                        y=df.Sample,
                        orientation="h",
                        x=df[column],
                        marker_color=color,
                    )
                    for name, column, color in zip(names, columns, colors)
                ]
            )
        else:

            # in multiqc_bowtie2, column paired_aligned_mate_one is multiplied by 2 as compared to
            # column PE one mate mapped uniquely. to be compatible with previous version,
            # we multiply by 2.
            columns = [
                "paired_aligned_one",  # mapped uniquely",
                "paired_aligned_discord_one",  # mapped discordantly uniquely
                "paired_aligned_mate_one",  # one mate mapped uniquely",
                "paired_aligned_multi",  # multimapped",
                "paired_aligned_mate_multi_halved",  # one mate multimapped
                "paired_aligned_mate_none_halved",  # neither mate aligned
            ]
            names = [
                "Mapped uniquely",
                "discordand unique mapping",
                "One mate multi mapped",
                "Multi mapped paired",
                "One mate multi-mapped",
                "Unaligned",
            ]
            self.df["paired_aligned_mate_one"] /= 2

            colors = ["#20568f", "#5c94ca", "#95ceff", "#f7a35c", "#ffeb75", "#981919"]

            S = self.df[columns].sum(axis=1).values
            df = self.df[columns].divide(S, axis=0) * 100
            df["Sample"] = self.df["Sample"]

            fig = go.Figure(
                data=[
                    go.Bar(
                        name=name,
                        y=df.Sample,
                        orientation="h",
                        x=df[column],
                        marker_color=color,
                    )
                    for name, column, color in zip(names, columns, colors)
                ],
                layout_xaxis_range=[0, 100],
            )

        fig.update_layout(barmode="stack", title="", xaxis_title="Mapping rate (percentage)")

        if html_code:
            return fig
        else:  # pragma: no cover
            fig.show()


class STAR(Reader):
    def __init__(self, filename):
        super().__init__(filename)

    def plot(self, html_code=False):
        import plotly.graph_objects as go

        fig = go.Figure()

        # get percentage

        df = self.df[["Sample"] + [x for x in self.df.columns if "_percent" in x]]

        fig = go.Figure(
            data=[
                go.Bar(
                    name="Uniquely mapped",
                    y=df.Sample,
                    orientation="h",
                    x=df["uniquely_mapped_percent"],
                    marker_color="#437bb1",
                ),
                go.Bar(
                    name="Mapped to multiple loci",
                    y=df.Sample,
                    orientation="h",
                    x=df["multimapped_percent"],
                    marker_color="#7cb5ec",
                ),
                go.Bar(
                    name="Mapped to too many loci",
                    y=df.Sample,
                    orientation="h",
                    x=df["multimapped_toomany_percent"],
                    marker_color="#f7a35c",
                ),
                go.Bar(
                    name="Unmapped: to many mismatches",
                    y=df.Sample,
                    orientation="h",
                    x=df["unmapped_mismatches_percent"],
                    marker_color="#e63491",
                ),
                go.Bar(
                    name="Unmapped: too short",
                    y=df.Sample,
                    orientation="h",
                    x=df["unmapped_tooshort_percent"],
                    marker_color="#b1084c",
                ),
                go.Bar(
                    name="Unmapped: other",
                    y=df.Sample,
                    orientation="h",
                    x=df["unmapped_other_percent"],
                    marker_color="#028ce2",
                ),
            ]
        )

        fig.update_layout(barmode="stack", title="Alignment Scores", xaxis_title="Mapping rate (percentage)")

        if html_code:
            return fig
        else:  # pragma: no cover
            fig.show()


class FeatureCounts(Reader):
    def __init__(self, filename):
        super().__init__(filename)

    def plot(self, html_code=False):
        import plotly.graph_objects as go

        fig = go.Figure()

        # feature counts create a set of key/values. sometimes there are no multimapped
        # and all are unmapped. e.g. B18355 and sometimes, unmapped is zero but lots
        # of multimapped e.g. B16162. The multiqc data set will store only a subset
        # and that may differ but coherent (the empty column is not there. So we end
        # up with two possiblities:

        # looks like with bowtie2, multi-mapped is zero, while with stars, unmapped is zero.
        # multiqc reports one or the other. In version <=0.17.2 with used mqc_feature_counts files
        # where one or the other is stored. in version >0.17.2 ze use multiqc_featurecounts.txt
        # where both are provided...so no way to distingiush them. We wil plot the two entries.

        columns = [
            "Assigned",
            "Unassigned_Unmapped",
            "Unassigned_MultiMapping",
            "Unassigned_NoFeatures",
            "Unassigned_Ambiguity",
        ]
        S = self.df[columns].sum(axis=1).values

        df = self.df[columns].divide(S, axis=0) * 100
        df["Sample"] = self.df["Sample"]

        colors = {
            "Assigned": "#7cb5ec",
            "Unassigned_Unmapped": "#434348",
            "Unassigned_MultiMapping": "#434348",
            "Unassigned_NoFeatures": "#90ed7d",
            "Unassigned_Ambiguity": "#f7a35c",
        }

        fig = go.Figure(
            data=[
                go.Bar(
                    name=name,
                    y=df.Sample,
                    orientation="h",
                    x=df[name],
                    marker_color=colors[name],
                )
                for name in columns
                if name in df.columns
            ]
        )

        fig.update_layout(
            barmode="stack", title="FeatureCounts: Assignments", xaxis_title="Annotation rate (percentage)"
        )

        if html_code:
            return fig
        else:  # pragma: no cover
            fig.show()
