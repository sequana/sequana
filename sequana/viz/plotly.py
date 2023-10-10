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


class BinaryPercentage:
    """
    Expects a dataframe with 2 columns. Their names are used for the labels.
    Indices of the dataframe is the sample name.

        from sequana.viz.plotly import HorizontalBar_BinaryPercentage
        hb = HorizontalBar_BinaryPercentage()
        hb.df = df
        hb.plot()
    """

    def __init__(self):
        self.df = None

    def plot_horizontal_bar(self, html_code=False, colors=["#389f1f", "#120946"]):
        import plotly.graph_objects as go

        fig = go.Figure()

        Xname = self.df.columns[0]
        Yname = self.df.columns[1]

        fig.add_trace(go.Bar(y=self.df.index, x=self.df[Yname], name=Yname, marker_color=colors[0], orientation="h"))
        fig.add_trace(go.Bar(y=self.df.index, x=self.df[Xname], name=Xname, marker_color=colors[1], orientation="h"))

        # Here we modify the tick angle of the xaxis, resulting in rotated labels.
        fig.update_layout(barmode="stack", xaxis_tickangle=-45, height=400)
        #    title="Mapped reads")

        if html_code:
            return fig
        else:
            fig.show()
