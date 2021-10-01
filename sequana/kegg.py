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
from sequana import sequana_data

from sequana.lazy import pandas as pd

import colorlog
logger = colorlog.getLogger(__name__)


__all__ = ["KEGGHelper"]


class KEGGHelper():
    """A simple class to build kegg information"""

    def __init__(self):


        self.df = pd.read_csv(sequana_data("kegg.csv"), index_col=0)
        self.df.fillna("", inplace=True)

    def build_csv(self, filename=None, Nmax=None):
        """rebuild the entire dataframe (1hour) and stores as attribute

        :param Nmax: for testing
        """
        logger.info("Retrieving the kegg organisms and their definitions")
        from bioservices import KEGG
        k = KEGG()
        results = []
        definition = []
        for i, item in enumerate(k.organismIds):
            results.append(k.parse(k.get(f"gn:{item}"))['NAME'])
            definition.append(k.parse(k.get(f"gn:{item}"))['ORG_CODE'])
            print(i, Nmax)
            if Nmax and i+1 >= Nmax:
                break

        results = [x[0] for x in results]
        IDs = [x.split(",")[0] for x in results]
        taxon = [x.split(",")[-1] for x in results]
        names = [x.split(",")[1].strip() if len(x.split(","))==3 else None for x in results]

        df = pd.DataFrame({'ID': IDs, 'taxon': taxon, 'name': names, 'def':definition })
        df = df.fillna("")
        df.columns = ['ID', 'taxon', 'shortname', 'definition']
        df['definition'] = [x.lower() for x in df.definition]
        df['shortname'] = [x.lower() for x in df.shortname]

        self.df = df
        if filename:
            df.to_csv(filename)

    def search(self, pattern):
        # if pattern is a string
        pattern = str(pattern)
        f1 = self.df[[True if pattern in x else False for x in self.df.definition]]
        f2 = self.df[[True if pattern in x else False for x in self.df.shortname]]
        f3 = self.df[[True if pattern in x else False for x in self.df.ID]]
        indices = list(f1.index) + list(f2.index) + list(f3.index)

        if len(indices) == 0:
            # maybe it is a taxon ID ?
            f4 = self.df[[True if pattern in str(x) else False for x in self.df.taxon]]
            indices = list(f4.index)

        results = self.df.loc[indices]
        return results

