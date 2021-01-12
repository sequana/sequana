# -*- coding: utf-8 -*-
#
#  This file is part of Sequana software
#
#  Copyright (c) 2016-2021 - Sequana Development Team
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
import uuid

from sequana.utils.datatables_js import DataTable


def df2html(df, name=None, dom='Brt'):
    """Simple wrapper to create HTML from dataframe"""

    if name is None:
        name = uuid.uuid1().time_low
        # looks like datatable does not like ID made of numbers, even in string
        # so we convert to ABCDEFGH values
        name = "".join([chr(65+int(x)) for x in str(name)])

    datatable = DataTable(df, name)
    datatable.datatable.datatable_options = {
            'pageLength': 15,
            'scrollCollapse': 'false',
            'dom': dom,
            'buttons': ['copy', 'csv']}
    js = datatable.create_javascript_function()
    html = datatable.create_datatable(float_format='%.6g')
    return js + html

