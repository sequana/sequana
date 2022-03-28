#
#  This file is part of Sequana software
#
#  Copyright (c) 2016-2022 - Sequana Development Team
#
#  Distributed under the terms of the 3-clause BSD license.
#  The full license is in the LICENSE file, distributed with this software.
#
#  website: https://github.com/sequana/sequana
#  documentation: http://sequana.readthedocs.io
#
##############################################################################
"""Tools to product json files for canvaJS plot.

"""
from collections import Counter
import json

from sequana.lazy import numpy as np

import colorlog

logger = colorlog.getLogger(__name__)


def list_to_json_for_barplot(l, logy=False):
    """Take a list and convert as json format for barplot with canvaJS.
    Return a string.

    :param: list l: list of int
    :param: bool logy: for logscale plot
    """
    counter = Counter(l)
    if logy:
        # tolist() necessary to convert numpy.float64 as float
        count = [{"x": int(key), "y": np.log10(value), "c": value} for key, value in counter.items()]
    else:
        count = [{"x": int(key), "y": value} for key, value in counter.items()]
    return json.dumps(count)
