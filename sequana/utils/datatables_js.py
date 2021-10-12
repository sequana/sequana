# coding: utf-8
#
#  This file is part of Sequana software
#
#  Copyright (c) 2016 - Sequana Development Team
#
#  File author(s):
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
""" Utilities to create a Jquery DataTable for your HTML file. 


.. autosummary::

    DataTableFunction
    DataTable


"""
from collections import OrderedDict

import colorlog
logger = colorlog.getLogger(__name__)



class DataTableFunction(object):
    """ Class that contains Jquery DataTables function and options.

    Example:

    ::

        import pandas as pd
        from sequana.utils import DataTableFunction

        df = pandas.read_csv('data.csv')
        datatable_js = DataTableFunction(df, 'data')
        datatable_js.datatable_options = {'pageLength': 15,
                                          'dom': 'Bfrtip',
                                          'buttons': ['copy', 'csv']}
        js = datatable_js.create_javascript_function()
        html_datatables = [DataTable(df, "data_{0}".format(i), datatable_js)
                           for i, df in enumerate(df_list)]

    Here, the datatable_options dictionary is used to fine tune the appearance
    of the table. 

    .. note::  DataTables add a number of elements around the table to control
        the table or show additional information about it. There are controlled
        by the order in the document (**DOM**) defined as a string made of 
        letters, each of them having a precise meaning. The order of the letter 
        is important. For instance if **B** is first, the buttons are put before 
        the table. If **B** is at the end, it is shown below the table. 
        Here are some of the valid letters and their meaning:

        - **B**: add the Buttons (copy/csv)
        - **i**: add *showing 1 to N of M entries*
        - **f**: add a search bar (**f** filtering)
        - **r**: processing display element
        - **t**: the table itself
        - **p**: pagination control

        Each option can be specified multiple times (with the exception of the 
        table itself).

    .. note:: other useful options are:

        - pageLength: 15
        - scrollX: "true"
        - paging: 15
        - buttons: ['copy', 'csv']

        Note that buttons can also be excel, pdf, print, ...



    All options of datatable:
        https://datatables.net/reference/option/
    """
    def __init__(self, df, html_id, index=False):
        """.. rubric:: contructor

        :param df: data frame.
        :param str html_id: the ID used in the HTML file.
        """
        self.index = index
        self._html_id = html_id
        self._datatable_options = dict()
        self._datatable_columns = self._set_datatable_columns(df)

    @property
    def html_id(self):
        """ Get the html_id, which cannot be set by the user after the
        instanciation of the class.
        """
        return self._html_id

    @property
    def datatable_options(self):
        """ Get, set or delete the DataTable options. Setter takes a dict as
        parameter with the desired options and updates the current dictionary.

        Example::

            datatable = DataTableFunction("tab")
            datatable.datatable_options = {'dom': 'Bfrtip',
                                           'buttons': ['copy', 'csv']}

        source: https://datatables.net/reference/option/
        """
        return self._datatable_options

    @datatable_options.setter
    def datatable_options(self, d):
        try:
            d['buttons'] = self._add_export_visible(d['buttons'])
        except KeyError:
            pass
        self._datatable_options.update(d)

    def _add_export_visible(self, buttons):
        """ Add option to disable the exporting of hidden columns
        """
        try:
            for b in buttons:
                b.update({'exportOptions': {'columns': ':visible'}})
        except AttributeError:
            buttons = [{'extend': b, 'exportOptions': {'columns': ':visible'}}
                      for b in buttons]
        return buttons

    @datatable_options.deleter
    def datatable_options(self):
        self._datatable_options = dict()

    @property
    def datatable_columns(self):
        """ Get datatable_columns dictionary. It is automatically set from the
        dataframe you want to plot.
        """
        return self._datatable_columns

    def _set_datatable_columns(self, df):
        """ Fill :attr:`DataTableFunction.datatable_columns` with header of
        :param:`DataTableFunction.df`.
        """
        from pandas import Series
        if isinstance(df, Series):
            return {}

        if self.index is True:
            columns = [""] + list(df.columns)
        else:
            columns = list(df.columns)
        column_dict = OrderedDict((name, dict()) for name in columns)
        return column_dict

    def create_javascript_function(self):
        """ Return javascript to create the DataTable.
        """
        js_function = """
<script type="text/javascript">
    function parseCsv_{0}(csv, id) {{
        Papa.parse(csv, {{
            comments: '#',
            delimiter: ',',
            header: true,
            dynamicTyping: true,
            error: function(reason) {{
                console.log(reason);
            }},
            complete: function(results) {{
                {1}
            }}
        }});
    }};
</script>
"""
        return js_function.format(self.html_id,
                                  self._create_datatable_option())

    def _create_datatable_option(self):
        """ Return DataTable options.
        """
        self.datatable_options['columns'] = self._create_columns_option()
        js = self._dict_to_string(self.datatable_options)
        js = "$(id).DataTable({{{0},data: results.data}});".format(js)
        return js

    def _create_columns_option(self):
        """ Return string well formated with all columns options.
        """
        js = [self._coloption_2_str(key, value) for key, value in
              self.datatable_columns.items()]
        return '[{0}]'.format(',\n'.join(js))

    def _coloption_2_str(self, name, options):
        s = "data:'{0}'".format(name)
        if options:
            s = "{0},\n{1}".format(s, self._dict_to_string(options))
        return '{{{0}}}'.format(s)

    def _dict_to_string(self, d):
        """ Convert dict to string for CanvasJS.

        Example:

        ::
            dico = {'key1': value1, 'key2': value2, 'key3': value3}
            print(CanvasJS._dict_to_string(dico))

            "key1:value1,key2:value2,key3:value3"
        """
        s = ['{0}:{1}'.format(key, self._check_type(value)) for key, value in
             d.items()]
        return ',\n'.join(s)

    def _check_type(self, value):
        """ Check value type to fill javascript sections. String must be
        surrounded by quotes and not boolean or integer.

        Javascript variable must not be surrounded by quotes. Custom variables
        start with 'data_'.
        """
        try:
            if not value.startswith(('true', 'false', 'function', '{', '[')):
                return "'{0}'".format(value)
        except AttributeError:
            return value
        return value

    def set_links_to_column(self, link_col, target_col, new_page=True):
        """Hide a column with urls and connect it with a column.

        :param str link_col: column with your URLs.
        :param str target_col: column to connect.
        """
        # hide the link column
        try:
            self.datatable_columns[link_col]['visible'] = 'false'
        except KeyError:

            keys = self.datatable_columns.keys() 
            logger.warning(f"KeyError: Column name '{target_col}' does not exist. Use one of {keys}")
        # function to add link
        if new_page is True:
            fct = """function(data, type, row, meta){{
                return '<a href="'+row.{0}+'" target="_blank">'+data+'</a>';
            }}
            """.format(link_col)
        else:
            fct = """function(data, type, row, meta){{
                return '<a href="'+row.{0}+'">'+data+'</a>';
            }}
            """.format(link_col)
        try:
            self.datatable_columns[target_col]['render'] = fct
        except KeyError:
            logger.warning("KeyError: Column name '{0}' does not exist."
                           .format(target_col))
            pass

    def set_tooltips_to_column(self, tooltips_col, target_col):
        """Hide a column with tooltips and connect it with a column.

        :param str tooltips_col: column with your tooltips.
        :param str target_col: column to connect.
        """
        # hide tooltips
        try:
            self.datatable_columns[tooltips_col]['visible'] = 'false'
        except KeyError:
            logger.warning("KeyError: Column name '{0}' does not exist."
                           .format(target_col))
            pass
        # function to add tooltips
        fct = """function(data, type, row, meta){{
            return '<a href="#" data-toggle="tooltip" title="'+row.{0}+'">'+data+'</a>';
        }}
        """.format(tooltips_col)
        try:
            self.datatable_columns[target_col]['render'] = fct
        except KeyError:
            logger.warning("KeyError: Column name '{0}' does not exist."
                           .format(target_col))
            pass


class DataTable(object):
    """ Class that contains html table which used a javascript function. 

    You  must add in your HTML file the JS function
    (:meth:`DataTable.create_javascript_function`) and the HTML code
    (:meth:`DataTable.create_datatable`).

    Example:

    ::

        df = pandas.read_csv('data.csv')
        datatable = DataTable(df, 'data')
        datatable.datatable.datatable_options = {'pageLength': 15,
                                                 'dom': 'Bfrtip',
                                                 'buttons': ['copy', 'csv']}
        js = datatable.create_javascript_function()
        html = datatable.create_datatable()

        # Second CSV file with same format
        df2 = pandas.read_csv('data2.csv')
        datatable2 = DataTable(df2, 'data2', datatable.datatable)
        html2 = datatable.create_datatable()

    The reason to include the JS manually is that you may include many HTML
    table but need to include the JS only once.

    """
    def __init__(self, df, html_id, datatable=None, index=False):
        """.. rubric:: contructor

        :param df: data frame.
        :param str html_id: the unique ID used in the HTML file.
        :param DataTableFunction datatable: javascript function to create the
            Jquery Datatables. If None, a :class:`DataTableFunction` is
            generated from the df.
        :param bool index: indicates whether the index dataframe should 
            be included in the CSV table
        """
        self.index = index
        self._df = df
        self._html_id = html_id
        if datatable:
            self.datatable = datatable
        else:
            self.datatable = DataTableFunction(df, html_id, index=index)

    def __len__(self):
        return len(self.df)

    @property
    def df(self):
        return self._df

    @property
    def html_id(self):
        return self._html_id

    def create_datatable(self, style="width:100%", **kwargs):
        """ Return string well formated to include in a HTML page.

        :param str style: CSS option of your table.
        :param dict kwargs: parameters of :meth:`pandas.DataFrame.to_csv`.
        """

        html = """
<script type="text/javascript">
    $(document).ready(function() {{
        var {0} = document.getElementById('csv_{0}').innerText;
        parseCsv_{1}({0}, '#table_{0}');
        {0} = null;
    }});
</script>
        """.format(self.html_id, self.datatable.html_id)
        html += self._create_hidden_csv(**kwargs)
        html += self._create_html_table(style)
        return html

    def _create_hidden_csv(self, **kwargs):
        """ Return the HTML code and the CSV code for your hidden CSV section.

        :param **dict kwargs: parameters of :meth:`pandas.DataFrame.to_csv`.
        """
        csv = self._df.to_csv(index=self.index, **kwargs)
        html = '<pre id="csv_{0}">{1}</pre>'.format(self.html_id, csv.strip())
        css = '<style>#csv_{0}{{display:none}}</style>'.format(self.html_id)
        return '{0}\n{1}\n'.format(css, html)

    def _create_html_table(self, style):
        """ Just for set some option and header.

        :param str style: css option of your table.
        """
        # set table id
        if style:
            style = 'style="{0}"'.format(style)
        html_table = (
            '<table id="table_{0}" class="display table text-center" {1}>'
            .format(self.html_id, style)
        )
        # create table's header
        th = '<th>{0}</th>'
        if self.index is True:
            header = [th.format("")]
            header += [th.format(name) for name in self.df]
        else:
            header = [th.format(name) for name in self.df]
        header = '<thead><tr>{0}</tr></thead>'.format("\n".join(header))
        html_table = """
    {0}
        {1}
    </table>
        """.format(html_table, header)
        return html_table

    def create_javascript_function(self):
        """ Generate the javascript function to create the DataTable in a HTML
        page.
        """
        return self.datatable.create_javascript_function()
