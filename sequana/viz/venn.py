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
from matplotlib_venn import venn2_unweighted, venn3_unweighted
from sequana.lazy import pylab


__all__ = ["plot_venn"]


def plot_venn(subsets, labels=None, title=None, ax=None, 
    alpha=0.8, weighted=False, colors=('r', 'b', 'y')):

    """ Plot venn diagramm according to number of groups.

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

    #pylab.clf()
    if len(subsets) == 2:
        from matplotlib_venn import venn2_unweighted, venn2_circles, venn2
        if weighted:
            venn_function = venn2
        else:
            venn_function = venn2_unweighted
        venn_circles = venn2_circles
    elif len(subsets) == 3:
        from matplotlib_venn import venn3_unweighted, venn3_circles, venn3
        if weighted:
            venn_function = venn3
        else:
            venn_function = venn3_unweighted
        venn_circles = venn3_circles
    else:
        raise IOError("Venn diagramm supports only 2 or 3 groups.")


    vf = venn_function(subsets, set_labels=labels, ax=ax, alpha=alpha,
        set_colors=colors)
    # works for weighted, nor for unweighted, so we draw the circles 
    # ourselfves
    if ax is None:
        ax = pylab.gca()

    for center, radius in zip(vf.centers, vf.radii):
        circle = pylab.Circle(center, radius=radius,
            linestyle="-", edgecolor="k", lw=2, facecolor='none', 
            alpha=1)
        ax.add_patch(circle)
    if title:
        pylab.title(title)
    return vf





























