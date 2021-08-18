# -*- coding: utf-8 -*-
#
#  This file is part of Sequana software
#
#  Copyright (c) 2016-2020 - Sequana Development Team
#
#  File author(s):
#      Thomas Cokelaer <thomas.cokelaer@pasteur.fr>
#      Etienne Kornobis <etienne.kornobis@pasteur.fr>
#
#  Distributed under the terms of the 3-clause BSD license.
#  The full license is in the LICENSE file, distributed with this software.
#
#  website: https://github.com/sequana/sequana
#  documentation: http://sequana.readthedocs.io
#
##############################################################################
"""Heatmap and dendograms"""
from sequana.lazy import pylab
from sequana.lazy import scipy
from sequana.lazy import numpy as np
from sequana.lazy import pandas as pd

# FIXME speed up this import
import scipy.cluster.hierarchy as hierarchy
import scipy.spatial.distance as distance

import easydev
from sequana.viz.linkage import Linkage
import matplotlib.patches as mpatches

__all__ = ["Heatmap", "Clustermap"]


try:
    import seaborn as sns
except ImportError:
    print("You should install seaborn to use sequana.viz.heatmap")


def get_heatmap_df():
    """a simple example to play with and perform test"""
    import pandas as pd

    df = pd.DataFrame(
        {
            "A": [1, 0, 1, 1],
            "B": [0.9, 0.1, 0.6, 1],
            "C": [0.5, 0.2, 0, 1],
            "D": [0.5, 0.2, 0, 1],
        }
    )
    return df


def get_clustermap_data():
    """a simple example to play with and perform test"""
    import pandas as pd

    df = pd.DataFrame(
        {
            "sample_1": [1, 0, 1, 1],
            "sample_2": [0.9, 0.1, 0.6, 1],
            "sample_3": [0.5, 0.2, 0, 1],
            "sample_4": [0.5, 0.2, 0, 1],
        },
        index=["gene_A", "gene_B", "gene_C", "gene_D"],
    )
    df_sample_groups = pd.DataFrame(
        {
            "sample_1": ["WT", "rep1"],
            "sample_2": ["WT", "rep2"],
            "sample_3": ["KO", "rep1"],
            "sample_4": ["KO", "rep2"],
        },
        index=["group1", "group2"],
    ).transpose()

    df_gene_groups = pd.DataFrame(
        {
            "gene_A": ["transcription", "GO1"],
            "gene_B": ["transcription", "GO2"],
            "gene_C": ["regulation", "GO1"],
            "gene_D": ["regulation", "GO2"],
        },
        index=["group1", "group2"],
    ).transpose()

    return df, df_sample_groups, df_gene_groups


# def heatmap(data, *args, **kargs):
#    """alias to Heatmap class"""
#    h = Heatmap(data, *args, **kargs)
#    h.plot()
#    return h


class Clustermap:
    """Heatmap and dendrograms based on seaborn Clustermap

    .. plot::
        :include-source:
        :width: 80%

        from sequana.viz.heatmap import Clustermap, get_clustermap_data
        df, sample_groups_df, gene_groups_df = get_clustermap_data()
        h = Clustermap(df, sample_groups_df=sample_groups_df, gene_groups_df=gene_groups_df)
        h.plot()

    """

    def __init__(
        self,
        data_df,
        sample_groups_df=None,
        sample_groups_sel=[],
        sample_groups_palette=sns.color_palette(desat=0.6),
        gene_groups_df=None,
        gene_groups_sel=[],
        gene_groups_palette=sns.color_palette(desat=0.6),
        yticklabels="auto",
        **kwargs
    ):  # annot):
        """.. rubric:: Constructor

        :param data_df: a dataframe.
        :param sample_groups_df: a dataframe with sample id as index (same as
            in data_df columns) and a group definition per column. Use to
            produce the x axis color groups.
        :param sample_group_sel: a list of the columns to select from the
            sample_groups_df.
        :param sample_groups_palette: the palette to use for sample color groups.
        :param gene_groups_df: a dataframe with gene id as index (same as in
            data_df columns) and a group definition per column. Use to produce
            the y axis color groups.
        :param gene_group_sel: a list of the columns to select from the gene_groups_df.
        :param gene_groups_palette: the palette to use for gene color groups.
        :param ytickslabels: "auto" for classical heatmap behaviour, [] for no
            ticks or a pandas Series giving the mapping between the index (gene
            names in data_df) and the gene names to be used for the heatmap
        :param kwargs: All other kwargs are passed to seaborn.Clustermap.

        """
        self.data_df = data_df
        self.sample_groups_df = sample_groups_df
        self.sample_groups_sel = sample_groups_sel
        self.sample_groups_palette = sample_groups_palette
        self.gene_groups_df = gene_groups_df
        self.gene_groups_sel = gene_groups_sel
        self.gene_groups_palette = gene_groups_palette
        self.kwargs = kwargs

        self.sample_groups_col_df, self.sample_color_dict = self._get_group_colors(
            sample_groups_df, sample_groups_sel, sample_groups_palette
        )
        self.gene_groups_col_df, self.gene_color_dict = self._get_group_colors(
            gene_groups_df, gene_groups_sel, gene_groups_palette
        )

        self.yticklabels = self._convert_gene_names(yticklabels)

    def _convert_gene_names(self, yticklabels):
        if isinstance(yticklabels, pd.Series):
            return yticklabels.loc[self.data_df.index]
        else:
            return yticklabels

    def _get_group_colors(self, groups_df, selection, palette):

        if isinstance(groups_df, pd.DataFrame):

            if selection:
                groups_df = groups_df.loc[:, selection]

            groups = [x for y in groups_df for x in groups_df[y].unique()]
            col_map = dict(zip(groups, palette))
            groups_col_df = groups_df.apply(lambda x: x.map(col_map))

            return groups_col_df, col_map
        else:
            return None, None

    def _do_legend(self, figure, color_dict, bbox_to_anchor):
        if color_dict:
            patches = [mpatches.Patch(color=c, label=l) for l, c in color_dict.items()]
            legend = pylab.legend(
                loc="upper center",
                handles=patches,
                bbox_to_anchor=bbox_to_anchor,
                frameon=True,
                title="Sample groups",
            )

            figure.add_artist(legend)

    def plot(self):
        cmap = sns.diverging_palette(220, 10, as_cmap=True)
        cmap.set_bad("grey", 1.0)

        p = sns.clustermap(
            self.data_df,
            cmap=cmap,
            col_colors=self.sample_groups_col_df,
            yticklabels=self.yticklabels,
            row_colors=self.gene_groups_col_df,
            **self.kwargs,
        )

        f = pylab.gca()

        self._do_legend(f, self.sample_color_dict, (12, 2))
        self._do_legend(f, self.gene_color_dict, (-2, -2))

        return p


class Heatmap(Linkage):
    """Heatmap and dendograms of an input matrix

    A heat map is an image representation of a matrix with a
    dendrogram added to the left side and to the top.  Typically,
    reordering of the rows and columns according to some set of values
    (row or column means) within the restrictions imposed by the
    dendrogram is carried out.


    .. plot::
        :include-source:
        :width: 80%

        from sequana.viz import heatmap
        df = heatmap.get_heatmap_df()
        h = heatmap.Heatmap(df)
        h.plot()


    side colors can be added:

        h = viz.Heatmap(df, col_side_colors=['r', 'g', 'b', 'y', 'k']); 
        h.category_column = category; 
        h.category_row = category

    where category is a dictionary with keys as df.columns and values as
    category defined by you. The number of colors in col_side_colors and
    row_side_colors should match the number of category 

    """

    def __init__(
        self,
        data=None,
        row_method="complete",
        column_method="complete",
        row_metric="euclidean",
        column_metric="euclidean",
        cmap="yellow_black_blue",
        col_side_colors=None,
        row_side_colors=None,
        verbose=True,
    ):
        """.. rubric:: constructor

        :param data: a dataframe or possibly a numpy matrix.
        :param row_method: complete by default
        :param column_method: complete by default. See linkage module for details
        :param row_metric: euclidean by default
        :param column_metric: euclidean by default
        :param cmap: colormap. any matplotlib accepted or combo of colors as
            defined in colormap package (pypi)
        :param col_side_colors:
        :param row_side_colors:


        """
        # should be a copy since it may be reshuffled ?
        try:
            if data is None and verbose is True:
                print("No data provided, please fill the `df` attribute manually")
            elif data is None:
                pass
            else:
                self._df = data.copy()
        except AttributeError as err:
            print("input must be a pandas data frame or numpy matrix")
            raise (err)

        self._row_method = row_method
        self._column_method = column_method

        self._column_metric = column_metric
        self._row_metric = row_metric

        # some default parameters
        self.cluster_criterion = "distance"
        self.params = easydev.AttrDict()
        self.params.col_side_colors = ["r", "g", "b", "y", "w", "k", "m"]
        self.params.row_side_colors = ["r", "g", "b", "y", "w", "k", "m"]
        self.params.cmap = cmap

        self.category_row = {}
        self.category_column = {}

        if col_side_colors:
            self.params.col_side_colors = col_side_colors
        if row_side_colors:
            self.params.row_side_colors = row_side_colors

    def _get_df(self):
        return self._df

    def _set_df(self, data):
        self._df = data.copy()

    df = property(_get_df, _set_df)
    frame = property(_get_df, _set_df)

    def _get_row_method(self):
        return self._row_method

    def _set_row_method(self, value):
        self.check_method(value)
        self._row_method = value

    row_method = property(_get_row_method, _set_row_method)

    def _get_col_method(self):
        return self._column_method

    def _set_col_method(self, value):
        self.check_method(value)
        self._column_method = value

    column_method = property(_get_col_method, _set_col_method)

    def _get_col_metric(self):
        return self._column_metric

    def _set_col_metric(self, value):
        self.check_metric(value)
        self._column_metric = value

    column_metric = property(_get_col_metric, _set_col_metric)

    def _get_row_metric(self):
        return self._row_metric

    def _set_row_metric(self, value):
        self.check_metric(value)
        self._row_metric = value

    row_metric = property(_get_row_metric, _set_row_metric)

    def plot(
        self,
        num=1,
        cmap=None,
        colorbar=True,
        vmin=None,
        vmax=None,
        colorbar_position="right",
        gradient_span="None",
        figsize=(12, 8),
        fontsize=None,
    ):
        """

        Using as input::

            df = pd.DataFrame({'A':[1,0,1,1],
                               'B':[.9,0.1,.6,1],
                            'C':[.5,.2,0,1],
                            'D':[.5,.2,0,1]})

        we can plot the heatmap + dendogram as follows::

            h = Heatmap(df)
            h.plot(vmin=0, vmax=1.1)


        .. plot::
            :include-source:
            :width: 80%

            from sequana.viz import heatmap
            df = heatmap.get_heatmap_df()
            h = heatmap.Heatmap(df)
            h.category_column['A'] = 1
            h.category_column['C'] = 1
            h.category_column['D'] = 2
            h.category_column['B'] = 2
            h.plot()


        """
        # save all parameters in a dict
        layout = {}

        if cmap is None:
            cmap = self.params.cmap
        try:
            import colormap

            cmap = colormap.cmap_builder(cmap)
        except:
            pass

        # keep track of row and column names for later.
        row_header = self.frame.index
        column_header = self.frame.columns

        import matplotlib

        # FIXME something clever for the fontsize
        if len(row_header) > 100 or len(column_header) > 100:
            matplotlib.rcParams["font.size"] = 6
        if len(row_header) > 50 or len(column_header) > 50:
            matplotlib.rcParams["font.size"] = 7
        if len(row_header) > 30 or len(column_header) > 30:
            matplotlib.rcParams["font.size"] = 8
        else:
            matplotlib.rcParams["font.size"] = 12
        if fontsize:
            matplotlib.rcParams["font.size"] = fontsize

        # scaling min/max range
        self.gradient_span = gradient_span  #'only_max'
        # min_to_max, min_to_max_centered, only_max, only_min

        if self.gradient_span == "min_to_max_centered":
            vmax = self.frame.max().max()
            vmin = self.frame.min().min()
            vmax = max([vmax, abs(vmin)])
            vmin = vmax * -1
        if self.gradient_span == "only_max":
            vmin = 0
            vmax = self.frame.max().max()
        if self.gradient_span == "only_min":
            vmin = self.frame.min().min()
            vmax = 0
        norm = matplotlib.colors.Normalize(vmin, vmax)

        # Scale the figure window size #
        fig = pylab.figure(num=num, figsize=figsize)
        fig.clf()

        # LAYOUT --------------------------------------------------
        # ax1 (dendrogram 1) on the left of the heatmap
        [ax1_x, ax1_y, ax1_w, ax1_h] = [0.05, 0.22, 0.2, 0.6]
        width_between_ax1_axr = 0.004
        # distance between the top color bar axis and the matrix
        height_between_ax1_axc = 0.004
        # Sufficient size to show
        color_bar_w = 0.015

        # axr, placement of row side colorbar
        # second to last controls the width of the side color bar - 0.015 when showing
        [axr_x, axr_y, axr_w, axr_h] = [0.31, 0.1, color_bar_w, 0.6]
        axr_x = ax1_x + ax1_w + width_between_ax1_axr
        axr_y = ax1_y
        axr_h = ax1_h
        width_between_axr_axm = 0.004

        # axc, placement of column side colorbar #
        # last one controls the hight of the top color bar - 0.015 when showing
        [axc_x, axc_y, axc_w, axc_h] = [0.4, 0.63, 0.5, color_bar_w]
        axc_x = axr_x + axr_w + width_between_axr_axm
        axc_y = ax1_y + ax1_h + height_between_ax1_axc
        height_between_axc_ax2 = 0.004

        # axm, placement of heatmap for the data matrix # why larger than 1?
        [axm_x, axm_y, axm_w, axm_h] = [0.4, 0.9, 2.5, 0.5]
        axm_x = axr_x + axr_w + width_between_axr_axm
        axm_y = ax1_y
        axm_h = ax1_h
        axm_w = axc_w

        # ax2 (dendrogram 2), on the top of the heatmap #
        [ax2_x, ax2_y, ax2_w, ax2_h] = [0.3, 0.72, 0.6, 0.15]
        ax2_x = axr_x + axr_w + width_between_axr_axm
        ax2_y = ax1_y + ax1_h + height_between_ax1_axc + axc_h + height_between_axc_ax2
        ax2_w = axc_w

        # axcb - placement of the color legend #
        if colorbar_position == "top left":
            [axcb_x, axcb_y, axcb_w, axcb_h] = [0.07, 0.88, 0.18, 0.09]
        elif colorbar_position == "right":
            [axcb_x, axcb_y, axcb_w, axcb_h] = [0.85, 0.2, 0.08, 0.6]
        else:
            raise ValueError("'top left' or 'right' accepted for now")

        # COMPUTATION DENDOGRAM 1 -------------------------------------
        if self.column_method:
            Y = self.linkage(
                self.frame.transpose(), self.column_method, self.column_metric
            )
            ax2 = fig.add_axes([ax2_x, ax2_y, ax2_w, ax2_h], frame_on=True)

            #     p=30,    truncate_mode=None,    color_threshold=None,    get_leaves=True,
            # orientation='top    labels=None,    count_sort=False,    distance_sort=False,
            #     show_leaf_counts=True,    no_plot=False,    no_labels=False,    leaf_font_size=None,
            #     leaf_rotation=None,    leaf_label_func=None,    show_contracted=False,
            #     link_color_func=None,    ax=None,    above_threshold_color='b',            #

            # color_threshold=0 and above_threshold_color='k' colors all
            # dendogram into black
            Z = hierarchy.dendrogram(
                Y,
                color_threshold=0,
                above_threshold_color="k",
                distance_sort="descending",
            )
            ind2 = hierarchy.fcluster(Y, 0.7 * max(Y[:, 2]), self.cluster_criterion)

            ax2.set_xticks([])
            ax2.set_yticks([])
            # apply the clustering for the array-dendrograms to the actual matrix data
            idx2 = Z["leaves"]
            self.frame = self.frame.iloc[:, idx2]
            # reorder the flat cluster to match the order of the leaves the dendrogram
            ind2 = ind2[idx2]
            layout["dendogram2"] = ax2
        else:
            idx2 = range(self.frame.shape[1])

        # COMPUTATION DENDOGRAM 2 ---------------------------------
        if self.row_method:
            Y = self.linkage(self.frame, self.row_method, self.row_metric)

            ax1 = fig.add_axes([ax1_x, ax1_y, ax1_w, ax1_h], frame_on=True)
            Z = hierarchy.dendrogram(
                Y,
                orientation="right",
                color_threshold=0,
                above_threshold_color="k",
                distance_sort="descending",
            )
            ind1 = hierarchy.fcluster(Y, 0.7 * max(Y[:, 2]), self.cluster_criterion)

            ax1.set_xticks([])
            ax1.set_yticks([])
            # apply the clustering for the array-dendrograms to the actual matrix data
            idx1 = Z["leaves"]
            self.frame = self.frame.iloc[idx1, :]
            # reorder the flat cluster to match the order of the leaves the dendrogram
            ind1 = ind1[idx1]
            layout["dendogram1"] = ax1
        else:
            idx1 = range(self.frame.shape[0])

        # HEATMAP itself
        axm = fig.add_axes([axm_x, axm_y, axm_w, axm_h])
        axm.imshow(
            self.frame,
            aspect="auto",
            origin="lower",
            interpolation="None",
            cmap=cmap,
            norm=norm,
        )
        axm.set_xticks([])
        axm.set_yticks([])
        layout["heatmap"] = axm

        # TEXT
        new_row_header = []
        new_column_header = []
        for i in range(self.frame.shape[0]):
            axm.text(
                self.frame.shape[1] - 0.5,
                i,
                "  " + str(row_header[idx1[i]]),
                verticalalignment="center",
            )
            new_row_header.append(
                row_header[idx1[i]] if self.row_method else row_header[i]
            )

        for i in range(self.frame.shape[1]):
            axm.text(
                i,
                -0.9,
                " " + str(column_header[idx2[i]]),
                rotation=90,
                verticalalignment="top",
                horizontalalignment="center",
            )
            new_column_header.append(
                column_header[idx2[i]] if self.column_method else column_header[i]
            )

        # CATEGORY column ------------------------------
        if self.category_column:
            axc = fig.add_axes([axc_x, axc_y, axc_w, axc_h])

            category_col = [self.category_column[self.df.columns[i]] for i in idx2]

            dc = np.array(category_col, dtype=int)
            dc.shape = (1, len(ind2))
            cmap_c = matplotlib.colors.ListedColormap(self.params.col_side_colors)
            axc.matshow(dc, aspect="auto", origin="lower", cmap=cmap_c)
            axc.set_xticks([])
            axc.set_yticks([])
            layout["category_column"] = axc

        # CATEGORY row -------------------------------
        if self.category_row:
            axr = fig.add_axes([axr_x, axr_y, axr_w, axr_h])
            # self.category_row must be a dictionary with names as found in the columns
            # of the dataframe.

            category_row = [self.category_row[self.df.index[i]] for i in idx1]

            dr = np.array(category_row, dtype=int)
            dr.shape = (len(category_row), 1)
            cmap_r = matplotlib.colors.ListedColormap(self.params.col_side_colors)
            axr.matshow(dr, aspect="auto", origin="lower", cmap=cmap_r)
            axr.set_xticks([])
            axr.set_yticks([])
            layout["category_row"] = axr

        # COLORBAR ----------------------
        if colorbar == True:
            axcb = fig.add_axes([axcb_x, axcb_y, axcb_w, axcb_h], frame_on=False)
            if colorbar_position == "right":
                orientation = "vertical"
            else:
                orientation = "horizontal"
            cb = matplotlib.colorbar.ColorbarBase(
                ax=axcb, cmap=cmap, norm=norm, orientation=orientation
            )
            # axcb.set_title("whatever")
            # max_cb_ticks = 5
            # axcb.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(max_cb_ticks))
            layout["colorbar"] = cb
            layout["colorbar_scalablemap"] = axcb

        #   could be useful
        self.d = {"ordered": self.frame.copy(), "rorder": idx1, "corder": idx2}

        return layout
