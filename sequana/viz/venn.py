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


__all__ = ["plot_venn"]


def plot_venn(compa_list, labels, title=None, ax=None):
    """ Plot venn diagramm according to number of groups.

    .. plot::
        :include-source:

        from sequana.rnadiff import plot_venn
        plot_venn((100,200,10), ("a", 'b'))


    .. plot::
        :include-source:

        from sequana.rnadiff import plot_venn
        plot_venn((100,200,300, 10,20,30,5), ("a", 'b', 'c')) 


    """

    if len(compa_list) == 2:
        venn_function = venn2_unweighted
    elif len(compa_list) == 3:
        venn_function = venn3_unweighted
    else:
        raise IOError("Venn diagramm supports only 2 or 3 groups.")

    venn_function(compa_list, set_labels=labels, ax=ax)
    # ax.set_title = title


        pass
