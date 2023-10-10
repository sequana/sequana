#
#  This file is part of Sequana software
#
#  Copyright (c) 2016 - Sequana Development Team
#
#  File author(s):
#      Thomas Cokelaer <thomas.cokelaer@pasteur.fr>
#
#  Distributed under the terms of the 3-clause BSD license.
#  The full license is in the LICENSE file, distributed with this software.
#
#  website: https://github.com/sequana/sequana
#  documentation: http://sequana.readthedocs.io
#
##############################################################################
"""Simple utilities for pandas"""

from sequana.lazy import pandas as pd


class PandasReader:
    def __init__(self, filename, sep="\t", columns=None, **kwargs):
        try:
            self.df = pd.read_csv(filename, sep=sep, **kwargs)
        except pd.errors.EmptyDataError:
            self.df = pd.DataFrame(columns=columns)

        if columns:
            self.df.columns = columns

        # If there is a header, let us strip it from spaces
        try:
            self.df.columns = [x.strip() for x in self.df.columns]
        except:
            pass

        # let us strip strings from spaces
        for x in self.df.columns:
            try:
                self.df[x] = self.df[x].str.strip()
            except:
                pass
