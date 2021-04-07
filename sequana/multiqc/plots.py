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

class Reader():
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

        x = self.df['not_aligned_percentage']
        compx = 100 - x
    
        fig.add_trace(go.Bar(
            y =self.df['Sample'],
            x=compx,
            name='Aligned',
            marker_color='#389f1f', orientation="h"
        ))
        fig.add_trace(go.Bar(
            y=self.df.Sample,
            x=self.df['not_aligned_percentage'],
            name='Not Aligned',
            marker_color='#120946', orientation='h'
        ))

        # Here we modify the tickangle of the xaxis, resulting in rotated labels.
        fig.update_layout(barmode='stack', xaxis_tickangle=-45, height=400,
            title="Mapping on ribosomal/contaminant")

        if html_code:
            return fig
        else:
            fig.show()
