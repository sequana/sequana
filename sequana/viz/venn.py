# -*- coding: utf-8 -*-
#
#  This file is part of Sequana software
#
#  Copyright (c) 2016 - Sequana Development Team
#
#  File author(s):
#      Thomas Cokelaer <thomas.cokelaer@pasteur.fr>
#      Etienne Kornobis <etienne.kornobis@pasteur.fr>,
#
#  Distributed under the terms of the 3-clause BSD license.
#  The full license is in the LICENSE file, distributed with this software.
#
#  website: https://github.com/sequana/sequana
#  documentation: http://sequana.readthedocs.io
#
##############################################################################
from sequana.lazy import pylab

__all__ = ["plot_venn"]


import math

# coding: utf-8
from itertools import chain

from matplotlib import colors

from sequana.lazy import pylab

default_colors = [
    [92 / 255.0, 192 / 255.0, 98 / 255.0, 0.5],
    [90 / 255.0, 155 / 255.0, 212 / 255.0, 0.5],
    [246 / 255.0, 236 / 255.0, 86 / 255.0, 0.6],
    [241 / 255.0, 90 / 255.0, 96 / 255.0, 0.4],
    [255 / 255.0, 117 / 255.0, 0 / 255.0, 0.3],
    [82 / 255.0, 82 / 255.0, 190 / 255.0, 0.2],
]


def draw_shape(ax, shape, params, fillcolor):
    """
    General function to draw a shape (ellipse, polygon).
    """
    import matplotlib.patches as patches

    if shape == "ellipse":
        shape_obj = patches.Ellipse(**params, color=fillcolor, lw=3)
    elif shape == "polygon":
        shape_obj = patches.Polygon(params["xy"], closed=True, color=fillcolor)
    else:
        raise ValueError(f"Unsupported shape: {shape}")
    ax.add_patch(shape_obj)


def draw_text(ax, x, y, text, color="black", fontsize=14, ha="center", va="center"):
    """
    Draw text on the plot.
    """
    ax.text(x, y, text, fontsize=fontsize, color=color, ha=ha, va=va)


def get_labels(data, fill=["number"]):
    """
    Generate a dictionary of labels for the Venn diagram.
    """
    N = len(data)
    sets_data = [set(d) for d in data]
    union_all = set(chain(*data))

    set_collections = {}
    for n in range(1, 2**N):
        key = bin(n)[2:].zfill(N)
        included = [sets_data[i] for i, b in enumerate(key) if b == "1"]
        excluded = [sets_data[i] for i, b in enumerate(key) if b == "0"]
        intersection = set.intersection(*included) if included else union_all
        difference = intersection - set.union(*excluded) if excluded else intersection
        set_collections[key] = difference

    labels = {k: "" for k in set_collections}
    for k, values in set_collections.items():
        if "logic" in fill:
            labels[k] += f"{k}: "
        if "number" in fill:
            labels[k] += str(len(values))
        if "percent" in fill:
            labels[k] += f" ({100 * len(values) / len(union_all):.1f}%)"
    return labels


def _venn_plot(labels, names, layout, colors=None, figsize=(9, 9), dpi=96, fontsize=14):
    """
    General function to create a Venn diagram with configurable layout.
    """
    colors = colors or [default_colors[i] for i in range(len(layout["shapes"]))]
    fig, ax = pylab.subplots(figsize=figsize, dpi=dpi)
    ax.set_axis_off()
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)

    # Draw shapes and labels
    for i, (shape, params) in enumerate(layout["shapes"]):
        draw_shape(ax, shape, params, colors[i])

    for text_params in layout["texts"]:
        draw_text(ax, **text_params)

    # Add group names
    for i, name in enumerate(names):
        draw_text(ax, **layout["names"][i], text=name, color="k")

    leg = ax.legend(names, loc="center left", bbox_to_anchor=(1.0, 0.5), fancybox=True)
    leg.get_frame().set_alpha(0.5)

    return fig, ax


def venn2(labels, names=["A", "B"], **options):
    """
    Plot a 2-set Venn diagram.
    """
    layout = {
        "shapes": [
            ("ellipse", {"xy": (0.375, 0.5), "width": 0.6, "height": 0.6, "angle": 0}),
            ("ellipse", {"xy": (0.625, 0.5), "width": 0.6, "height": 0.6, "angle": 0}),
        ],
        "texts": [
            {"x": 0.74, "y": 0.50, "text": labels.get("01", "")},
            {"x": 0.26, "y": 0.50, "text": labels.get("10", "")},
            {"x": 0.50, "y": 0.50, "text": labels.get("11", "")},
        ],
        "names": [
            {"x": 0.20, "y": 0.76, "ha": "right", "va": "bottom"},
            {"x": 0.80, "y": 0.76, "ha": "left", "va": "bottom"},
        ],
    }
    return _venn_plot(labels, names, layout, **options)


def venn3(labels, names=["A", "B", "C"], **options):
    """
    Plot a 3-set Venn diagram.
    """
    layout = {
        "shapes": [
            ("ellipse", {"xy": (0.333, 0.633), "width": 0.5, "height": 0.5, "angle": 0}),
            ("ellipse", {"xy": (0.666, 0.633), "width": 0.5, "height": 0.5, "angle": 0}),
            ("ellipse", {"xy": (0.500, 0.310), "width": 0.5, "height": 0.5, "angle": 0}),
        ],
        "texts": [
            {"x": 0.50, "y": 0.27, "text": labels.get("001", "")},
            {"x": 0.73, "y": 0.65, "text": labels.get("010", "")},
            {"x": 0.61, "y": 0.46, "text": labels.get("011", "")},
            {"x": 0.27, "y": 0.65, "text": labels.get("100", "")},
            {"x": 0.39, "y": 0.46, "text": labels.get("101", "")},
            {"x": 0.50, "y": 0.65, "text": labels.get("110", "")},
            {"x": 0.50, "y": 0.51, "text": labels.get("111", "")},
        ],
        "names": [
            {"x": 0.15, "y": 0.87, "ha": "right", "va": "bottom"},
            {"x": 0.85, "y": 0.87, "ha": "left", "va": "bottom"},
            {"x": 0.50, "y": 0.02, "ha": "center", "va": "top"},
        ],
    }
    return _venn_plot(labels, names, layout, **options)


def venn4(labels, names=["A", "B", "C", "D"], **options):
    """
    Plot a 4-set Venn diagram.
    """
    layout = {
        "shapes": [
            ("ellipse", {"xy": (0.35, 0.4), "width": 0.72, "height": 0.45, "angle": 140}),
            ("ellipse", {"xy": (0.45, 0.5), "width": 0.72, "height": 0.45, "angle": 140}),
            ("ellipse", {"xy": (0.544, 0.5), "width": 0.72, "height": 0.45, "angle": 40}),
            ("ellipse", {"xy": (0.644, 0.4), "width": 0.72, "height": 0.45, "angle": 40}),
        ],
        "texts": [
            {"x": 0.85, "y": 0.42, "text": labels.get("0001", "")},
            {"x": 0.68, "y": 0.72, "text": labels.get("0010", "")},
            {"x": 0.77, "y": 0.59, "text": labels.get("0011", "")},
            {"x": 0.32, "y": 0.72, "text": labels.get("0100", "")},
            {"x": 0.71, "y": 0.30, "text": labels.get("0101", "")},
            {"x": 0.50, "y": 0.66, "text": labels.get("0110", "")},
            {"x": 0.65, "y": 0.50, "text": labels.get("0111", "")},
            {"x": 0.14, "y": 0.42, "text": labels.get("1000", "")},
            {"x": 0.50, "y": 0.17, "text": labels.get("1001", "")},
            {"x": 0.29, "y": 0.30, "text": labels.get("1010", "")},
            {"x": 0.39, "y": 0.24, "text": labels.get("1011", "")},
            {"x": 0.23, "y": 0.59, "text": labels.get("1100", "")},
            {"x": 0.61, "y": 0.24, "text": labels.get("1101", "")},
            {"x": 0.35, "y": 0.50, "text": labels.get("1110", "")},
            {"x": 0.50, "y": 0.38, "text": labels.get("1111", "")},
        ],
        "names": [
            {"x": 0.13, "y": 0.18, "ha": "center", "va": "center"},
            {"x": 0.18, "y": 0.83, "ha": "center", "va": "center"},
            {"x": 0.82, "y": 0.83, "ha": "center", "va": "center"},
            {"x": 0.87, "y": 0.18, "ha": "center", "va": "center"},
        ],
    }
    return _venn_plot(labels, names, layout, **options)


def venn5(labels, names=["A", "B", "C", "D", "E"], **options):
    """
    Plot a 5-set Venn diagram.
    """
    layout = {
        "shapes": [
            ("ellipse", {"xy": (0.428, 0.449), "width": 0.87, "height": 0.5, "angle": 155}),
            ("ellipse", {"xy": (0.469, 0.543), "width": 0.87, "height": 0.5, "angle": 82}),
            ("ellipse", {"xy": (0.558, 0.523), "width": 0.87, "height": 0.5, "angle": 10}),
            ("ellipse", {"xy": (0.578, 0.432), "width": 0.87, "height": 0.5, "angle": 118}),
            ("ellipse", {"xy": (0.489, 0.383), "width": 0.87, "height": 0.5, "angle": 46}),
        ],
        "texts": [
            {"x": 0.27, "y": 0.11, "text": labels.get("00001", "")},
            {"x": 0.72, "y": 0.11, "text": labels.get("00010", "")},
            {"x": 0.55, "y": 0.13, "text": labels.get("00011", "")},
            {"x": 0.91, "y": 0.58, "text": labels.get("00100", "")},
            {"x": 0.78, "y": 0.64, "text": labels.get("00101", "")},
            {"x": 0.84, "y": 0.41, "text": labels.get("00110", "")},
            {"x": 0.76, "y": 0.55, "text": labels.get("00111", "")},
            {"x": 0.51, "y": 0.90, "text": labels.get("01000", "")},
            {"x": 0.39, "y": 0.15, "text": labels.get("01001", "")},
            {"x": 0.42, "y": 0.78, "text": labels.get("01010", "")},
            {"x": 0.50, "y": 0.15, "text": labels.get("01011", "")},
            {"x": 0.67, "y": 0.76, "text": labels.get("01100", "")},
            {"x": 0.70, "y": 0.71, "text": labels.get("01101", "")},
            {"x": 0.51, "y": 0.74, "text": labels.get("01110", "")},
            {"x": 0.64, "y": 0.67, "text": labels.get("01111", "")},
            {"x": 0.10, "y": 0.61, "text": labels.get("10000", "")},
            {"x": 0.20, "y": 0.31, "text": labels.get("10001", "")},
            {"x": 0.76, "y": 0.25, "text": labels.get("10010", "")},
            {"x": 0.65, "y": 0.23, "text": labels.get("10011", "")},
            {"x": 0.18, "y": 0.50, "text": labels.get("10100", "")},
            {"x": 0.21, "y": 0.37, "text": labels.get("10101", "")},
            {"x": 0.81, "y": 0.37, "text": labels.get("10110", "")},
            {"x": 0.74, "y": 0.40, "text": labels.get("10111", "")},
            {"x": 0.27, "y": 0.70, "text": labels.get("11000", "")},
            {"x": 0.34, "y": 0.25, "text": labels.get("11001", "")},
            {"x": 0.33, "y": 0.72, "text": labels.get("11010", "")},
            {"x": 0.51, "y": 0.22, "text": labels.get("11011", "")},
            {"x": 0.25, "y": 0.58, "text": labels.get("11100", "")},
            {"x": 0.28, "y": 0.39, "text": labels.get("11101", "")},
            {"x": 0.36, "y": 0.66, "text": labels.get("11110", "")},
            {"x": 0.51, "y": 0.47, "text": labels.get("11111", "")},
        ],
        "names": [
            {"x": 0.0, "y": 0.62, "ha": "center", "va": "center"},
            {"x": 0.52, "y": 1.0, "ha": "center", "va": "center"},
            {"x": 0.99, "y": 0.74, "ha": "center", "va": "center"},
            {"x": 0.78, "y": 0.0, "ha": "center", "va": "center"},
            {"x": 0.22, "y": 0.0, "ha": "center", "va": "center"},
        ],
    }
    return _venn_plot(labels, names, layout, **options)


def venn6(labels, names=["A", "B", "C", "D"], **options):
    """
    Plot a 6-set Venn diagram.
    """

    layout = {
        "shapes": [
            ("polygon", {"xy": [(0.637, 0.921), (0.649, 0.274), (0.188, 0.667)]}),
            ("polygon", {"xy": [(0.981, 0.769), (0.335, 0.191), (0.393, 0.671)]}),
            ("polygon", {"xy": [(0.941, 0.397), (0.292, 0.475), (0.456, 0.747)]}),
            ("polygon", {"xy": [(0.662, 0.119), (0.316, 0.548), (0.662, 0.700)]}),
            ("polygon", {"xy": [(0.309, 0.081), (0.374, 0.718), (0.681, 0.488)]}),
            ("polygon", {"xy": [(0.016, 0.626), (0.726, 0.687), (0.522, 0.327)]}),
        ],
        "texts": [
            {"x": 0.212, "y": 0.562, "text": labels.get("000001", "")},
            {"x": 0.430, "y": 0.249, "text": labels.get("000010", "")},
            {"x": 0.356, "y": 0.444, "text": labels.get("000011", "")},
            {"x": 0.609, "y": 0.255, "text": labels.get("000100", "")},
            {"x": 0.323, "y": 0.546, "text": labels.get("000101", "")},
            {"x": 0.513, "y": 0.316, "text": labels.get("000110", "")},
            {"x": 0.523, "y": 0.348, "text": labels.get("000111", "")},
            {"x": 0.747, "y": 0.458, "text": labels.get("001000", "")},
            {"x": 0.325, "y": 0.492, "text": labels.get("001001", "")},
            {"x": 0.670, "y": 0.481, "text": labels.get("001010", "")},
            {"x": 0.359, "y": 0.478, "text": labels.get("001011", "")},
            {"x": 0.653, "y": 0.444, "text": labels.get("001100", "")},
            {"x": 0.344, "y": 0.526, "text": labels.get("001101", "")},
            {"x": 0.653, "y": 0.466, "text": labels.get("001110", "")},
            {"x": 0.363, "y": 0.503, "text": labels.get("001111", "")},
            {"x": 0.750, "y": 0.616, "text": labels.get("010000", "")},
            {"x": 0.682, "y": 0.654, "text": labels.get("010001", "")},
            {"x": 0.402, "y": 0.310, "text": labels.get("010010", "")},
            {"x": 0.392, "y": 0.421, "text": labels.get("010011", "")},
            {"x": 0.653, "y": 0.691, "text": labels.get("010100", "")},
            {"x": 0.651, "y": 0.644, "text": labels.get("010101", "")},
            {"x": 0.490, "y": 0.340, "text": labels.get("010110", "")},
            {"x": 0.468, "y": 0.399, "text": labels.get("010111", "")},
            {"x": 0.692, "y": 0.545, "text": labels.get("011000", "")},
            {"x": 0.666, "y": 0.592, "text": labels.get("011001", "")},
            {"x": 0.665, "y": 0.496, "text": labels.get("011010", "")},
            {"x": 0.374, "y": 0.470, "text": labels.get("011011", "")},
            {"x": 0.653, "y": 0.537, "text": labels.get("011100", "")},
            {"x": 0.652, "y": 0.579, "text": labels.get("011101", "")},
            {"x": 0.653, "y": 0.488, "text": labels.get("011110", "")},
            {"x": 0.389, "y": 0.486, "text": labels.get("011111", "")},
            {"x": 0.553, "y": 0.806, "text": labels.get("100000", "")},
            {"x": 0.313, "y": 0.604, "text": labels.get("100001", "")},
            {"x": 0.388, "y": 0.694, "text": labels.get("100010", "")},
            {"x": 0.375, "y": 0.633, "text": labels.get("100011", "")},
            {"x": 0.605, "y": 0.359, "text": labels.get("100100", "")},
            {"x": 0.334, "y": 0.555, "text": labels.get("100101", "")},
            {"x": 0.582, "y": 0.397, "text": labels.get("100110", "")},
            {"x": 0.542, "y": 0.372, "text": labels.get("100111", "")},
            {"x": 0.468, "y": 0.708, "text": labels.get("101000", "")},
            {"x": 0.355, "y": 0.572, "text": labels.get("101001", "")},
            {"x": 0.420, "y": 0.679, "text": labels.get("101010", "")},
            {"x": 0.375, "y": 0.597, "text": labels.get("101011", "")},
            {"x": 0.641, "y": 0.436, "text": labels.get("101100", "")},
            {"x": 0.348, "y": 0.538, "text": labels.get("101101", "")},
            {"x": 0.635, "y": 0.453, "text": labels.get("101110", "")},
            {"x": 0.370, "y": 0.548, "text": labels.get("101111", "")},
            {"x": 0.594, "y": 0.689, "text": labels.get("110000", "")},
            {"x": 0.579, "y": 0.670, "text": labels.get("110001", "")},
            {"x": 0.398, "y": 0.670, "text": labels.get("110010", "")},
            {"x": 0.395, "y": 0.653, "text": labels.get("110011", "")},
            {"x": 0.633, "y": 0.682, "text": labels.get("110100", "")},
            {"x": 0.616, "y": 0.656, "text": labels.get("110101", "")},
            {"x": 0.587, "y": 0.427, "text": labels.get("110110", "")},
            {"x": 0.526, "y": 0.415, "text": labels.get("110111", "")},
            {"x": 0.495, "y": 0.677, "text": labels.get("111000", "")},
            {"x": 0.505, "y": 0.648, "text": labels.get("111001", "")},
            {"x": 0.428, "y": 0.663, "text": labels.get("111010", "")},
            {"x": 0.430, "y": 0.631, "text": labels.get("111011", "")},
            {"x": 0.639, "y": 0.524, "text": labels.get("111100", "")},
            {"x": 0.591, "y": 0.604, "text": labels.get("111101", "")},
            {"x": 0.622, "y": 0.477, "text": labels.get("111110", "")},
            {"x": 0.501, "y": 0.523, "text": labels.get("111111", "")},
        ],
        "names": [
            {"x": 0.674, "y": 0.824, "ha": "center", "va": "center"},
            {"x": 0.747, "y": 0.751, "ha": "center", "va": "center"},
            {"x": 0.739, "y": 0.396, "ha": "center", "va": "center"},
            {"x": 0.700, "y": 0.247, "ha": "center", "va": "center"},
            {"x": 0.291, "y": 0.255, "ha": "center", "va": "center"},
            {"x": 0.203, "y": 0.484, "ha": "center", "va": "center"},
        ],
    }
    return _venn_plot(labels, names, layout, **options)


def plot_venn(subsets, labels=None, title=None, ax=None, alpha=0.8, weighted=False, colors=("r", "b", "y")):
    """Plot venn diagramm according to number of groups.

    :param subsets: This parameter may be (1) a dict, providing sizes of
        three diagram regions. The regions are identified via two-letter
        binary codes ('10', '01', and '11'), hence a valid set could look like:
        {'10': 10, '01': 20, '11': 40}. Unmentioned codes are considered to map to 0.
        (2) a list (or a tuple) with three numbers, denoting the sizes of the
        regions in the following order: (10, 01, 11) and (3) a list containing
        the subsets of values.

    The subsets can be a list (or a tuple) containing two set objects. For
    instance:

    .. plot::
        :include-source:

        from sequana.viz.venn import plot_venn
        A = set([1,2,3,4,5,6,7,8,9])
        B = set([            7,8,9,10,11])
        plot_venn((A, B), labels=("A", "B"))

    This is the unweighted version by default meaning all circles have the same
    size. If you prefer to have circle scaled to the
    size of the sets, add the relevant parameter as follows:

    .. plot::
        :include-source:

        from sequana.viz.venn import plot_venn
        A = set([1,2,3,4,5,6,7,8,9])
        B = set([            7,8,9,10,11])
        plot_venn((A, B), labels=("A", "B"), weighted=True)

    Similarly for 3 sets, a Venn diagram can be represented as follows. Note
    here that we also use the *title* parameter:

    .. plot::
        :include-source:

        from sequana.viz.venn import plot_venn

        A = set([1,2,3,4,5,6,7,8,9])
        B = set([      4,5,6,7,8,9,10,11,12,13])
        C = set([   3,4,5,6,7,8,9])
        plot_venn((A, B, C), labels=("A", "B", "C"), title="my Venn3 diagram")

    Input can be a list/tuple of 2 or 3 sets as described above.

    """
    if isinstance(subsets, (list, tuple)):
        n_sets = len(subsets)
        data = get_labels(subsets)
    elif isinstance(subsets, dict):
        n_sets = len(list(subsets.keys())[0])
        data = subsets
    venn_func_name = f"venn{n_sets}"

    try:
        venn_func = globals().get(venn_func_name)
        if not venn_func:
            raise ValueError(f"Venn diagram for {n_sets} sets is not supported.")
        venn = venn_func(data, names=labels)
        if title:
            pylab.title(title)
    except Exception as e:
        print(f"Error: {e}")
