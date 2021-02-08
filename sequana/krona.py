# -*- coding: utf-8 -*-
#
#  This file is part of Sequana software
#
#  Copyright (c) 2016 - Sequana Development Team
#
#  File author(s):
#      Thomas Cokelaer <thomas.cokelaer@pasteur.fr>
#      Dimitri Desvillechabrol <dimitri.desvillechabrol@pasteur.fr>,
#          <d.desvillechabrol@gmail.com>
#
#  Distributed under the terms of the 3-clause BSD license.
#  The full license is in the LICENSE file, distributed with this software.
#
#  website: https://github.com/sequana/sequana
#  documentation: http://sequana.readthedocs.io
#
##############################################################################
import collections

from sequana.lazy import pandas as pd
import colorlog
logger = colorlog.getLogger(__name__)




__all__ = ['KronaMerger']


class KronaMerger(collections.Counter):
    """Utility to merge two Krona files

    Imagine those two files (formatted for Krona; first column is a counter)::

        14011   Bacteria    Proteobacteria  species1
        591 Bacteria    Proteobacteria  species4

        184 Bacteria    Proteobacteria  species3
        132 Bacteria    Proteobacteria  species2
        32  Bacteria    Proteobacteria  species1

    You can merge the two files. The first and last lines correspond to the same
    taxon (species1) so we should end up with a new Krona file with 4 lines
    only.

    The test files are available within Sequana as test_krona_k1.tsv
    and test_krona_k2.tsv::

        from sequana import KronaMerger, sequana_data
        k1 = KronaMerger(sequana_data("test_krona_k1.tsv"))
        k2 = KronaMerger(sequana_data("test_krona_k2.tsv"))
        k1 += k2
        # Save the results. Note that it must be tabulated for Krona external usage
        k1.to_tsv("new.tsv")


    .. warning:: separator must be tabulars

    """
    def __init__(self, filename):
        """.. rubric:: constructor

        :param str filename:

        """
        super(KronaMerger, self).__init__()
        self.filename = filename
        self._read()

    def _read(self):
        with open(self.filename, "r") as fin:
            for line in fin.readlines():
                count, name = line.split("\t", 1)
                count = int(count)
                self[name] += count

    def to_tsv(self, output_filename):
        """Save the content into a new file in TSV format"""
        assert output_filename.endswith('.tsv')
        labels = []
        counts = []
        for k, count in self.items():
            labels.append(k)
            counts.append(count)
        df = pd.DataFrame({'label':labels, 'count':counts})
        df = df[['count', 'label']]
        df['label'] = df['label'].apply(lambda x: x.strip())
        try:
            df.sort_values("count", inplace=True, ascending=False)
        except:
            df.sort("count", inplace=True, ascending=False)
        df.to_csv(output_filename, sep="\t", index=None, header=None)
        return df
