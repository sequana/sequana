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


def df2html(df, name=None, dom='Brt', show_index=False, pageLength=15):
    """Simple wrapper to create HTML from dataframe

    If a columns ends in _links and a name_links exists, then the columns name 
    will be shown with the clickable name_links.
    """

    if name is None:
        name = uuid.uuid1().time_low
        # looks like datatable does not like ID made of numbers, even in string
        # so we convert to ABCDEFGH values
        name = "".join([chr(65+int(x)) for x in str(name)])

    datatable = DataTable(df, name, index=show_index)
    datatable.datatable.datatable_options = {
            'pageLength': pageLength,
            'scrollCollapse': 'false',
            'dom': dom,
            'buttons': ['copy', 'csv']}

    # identify links (columns ending in _links)
    for column in df.columns:
        if column.endswith('_links'):
            prefix = column.replace('_links', '')
            if prefix in df.columns:
                datatable.datatable.set_links_to_column(column, prefix)



    js = datatable.create_javascript_function()
    html = datatable.create_datatable(float_format='%.6g')
    return js + html

