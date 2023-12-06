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
        else:
            fig.show()


class Bowtie2(Reader):
    def __init__(self, filename):
        super().__init__(filename)

    def plot(self, html_code=False):
        import plotly.graph_objects as go

        fig = go.Figure()

        # get percentage instead of counts
        if "SE mapped uniquely" in self.df:

            S = self.df[["SE mapped uniquely", "SE multimapped", "SE not aligned"]].sum(axis=1).values
            df = self.df[["SE mapped uniquely", "SE multimapped", "SE not aligned"]].divide(S, axis=0) * 100
            df["Sample"] = self.df["Sample"]

            fig = go.Figure(
                data=[
                    go.Bar(
                        name="SE mapped uniquely",
                        y=df.Sample,
                        orientation="h",
                        x=df["SE mapped uniquely"],
                        marker_color="#17478f",
                    ),
                    go.Bar(
                        name="SE multimapped",
                        y=df.Sample,
                        orientation="h",
                        x=df["SE multimapped"],
                        marker_color="#e2780d",
                    ),
                    go.Bar(
                        name="SE not alinged",
                        y=df.Sample,
                        orientation="h",
                        x=df["SE not aligned"],
                        marker_color="#9f1416",
                    ),
                ]
            )
        else:
            columns = [
                "PE mapped uniquely",
                "PE mapped discordantly uniquely",
                "PE one mate mapped uniquely",
                "PE multimapped",
                "PE one mate multimapped",
                "PE neither mate aligned",
            ]
            colors = ["#20568f", "#5c94ca", "#95ceff", "#f7a35c", "#ffeb75", "#981919"]

            S = self.df[columns].sum(axis=1).values
            df = self.df[columns].divide(S, axis=0) * 100
            df["Sample"] = self.df["Sample"]

            fig = go.Figure(
                data=[
                    go.Bar(
                        name=column,
                        y=df.Sample,
                        orientation="h",
                        x=df[column],
                        marker_color=color,
                    )
                    for column, color in zip(columns, colors)
                ],
                layout_xaxis_range=[0, 100],
            )

        fig.update_layout(barmode="stack", title="", xaxis_title="Mapping rate (percentage)")

        if html_code:
            return fig
        else:
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
                    name="uniquely mapped",
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
        else:
            fig.show()


class FeatureCounts(Reader):
    def __init__(self, filename):
        super().__init__(filename)

    def plot(self, html_code=False):
        import plotly.graph_objects as go

        fig = go.Figure()

        columns = ["Assigned", "Unassigned: Unmapped", "Unassigned: No Features", "Unassigned: Ambiguity"]
        S = self.df[columns].sum(axis=1).values
        df = self.df[columns].divide(S, axis=0) * 100
        df["Sample"] = self.df["Sample"]

        fig = go.Figure(
            data=[
                go.Bar(
                    name="Assigned",
                    y=df.Sample,
                    orientation="h",
                    x=df["Assigned"],
                    marker_color="#7cb5ec",
                ),
                go.Bar(
                    name="Unassigned: Unmapped",
                    y=df.Sample,
                    orientation="h",
                    x=df["Unassigned: Unmapped"],
                    marker_color="#434348",
                ),
                go.Bar(
                    name="Unassigned: No Features",
                    y=df.Sample,
                    orientation="h",
                    x=df["Unassigned: No Features"],
                    marker_color="#90ed7d",
                ),
                go.Bar(
                    name="Unassigned: Ambiguity",
                    y=df.Sample,
                    orientation="h",
                    x=df["Unassigned: Ambiguity"],
                    marker_color="#f7a35c",
                ),
            ]
        )

        fig.update_layout(
            barmode="stack", title="FeatureCounts: Assignments", xaxis_title="Annotation rate (percentage)"
        )

        if html_code:
            return fig
        else:
            fig.show()
