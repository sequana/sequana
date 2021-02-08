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
"""Sequana class to plot a CanvasJS linegraph from an embedded csv file.
"""
from sequana.plots.canvasjs_base import CanvasJS

import colorlog
logger = colorlog.getLogger(__name__)



class CanvasJSLineGraph(CanvasJS):
    """ Class to create a CanvasJS linegraphe for an HTML page. It creates a
    hidden pre section with your CSV. It is necessary embedded because browsers
    forbid the reading of data stored locally. Your html page need CanvasJS and
    PapaParse.
    """
    def __init__(self, csv, html_id, x_column, y_columns):
        """.. rubric:: constructor

        :param str csv: data as CSV format.
        :param str html_id: the ID used in your html. All function
            will have this tag.
        :param str x_column: column used as axis X for plot.
        :param list y_columns: colums used as axis Y for plot.
        """
        super().__init__(html_id)
        self.csv = csv.strip()
        self.html_id = html_id
        self.x_column = x_column
        self.y_columns = y_columns
        # create hidden csv
        self.html_cjs = self._create_hidden_csv()
        self.html_cjs += '<script type="text/javascript">{0}\n'.format(
            self._create_js_csv_parser())

    def _create_hidden_csv(self):
        """ Return the HTML code and the CSV code for your hidden CSV section.
        """
        html = '<pre id="{0}">{1}</pre>'.format(self.html_id, self.csv)
        css = '<style>#{0}{{display:none}}</style>'.format(self.html_id)
        return '{0}\n{1}\n'.format(html, css)

    def _create_js_csv_parser(self):
        """ Create the efficient javascript csv parser with PapaParse.
        """
        # Create variable name
        variable = ["data_{0}".format(name) for name in self.y_columns]
        # Create variable for javascript array
        init_var = ["var {0} = [];".format(name) for name in variable]
        init_var = "\n".join(init_var)
        # Fill arrays
        fill_array = ["{0}.push({{x: curow.{1}, y: curow.{2}}});".format(v,
                      self.x_column, n) for n, v in
                      zip(self.y_columns, variable)]
        fill_array = "\n".join(fill_array)
        self.variables = ", ".join(variable)
        function = """
    function processData_{0}(csv) {{
        {1}
        var curow;
        Papa.parse(csv, {{
            comments: '#',
            delimiter: ',',
            header: true,
            dynamicTyping: true,
            error: function(reason) {{
                console.log(reason);
            }},
            step: function(row){{
                curow = row.data[0];
                {2}
            }},
            complete: function(results) {{
                drawChart_{0}({3});
            }}
        }});
    }};
        """.format(self.html_id, init_var, fill_array, self.variables)
        self.data_section = [{'dataPoints': var} for var in variable]
        return function

    def set_axis_x(self, axis_attr=dict()):
        """ Method to configure X axis of the line graph.

        :param dict axis_attr: dictionary with canvasjs axisX
            Attributes.

        axisX: http://canvasjs.com/docs/charts/chart-options/axisx/

        Example:

        ::

            line_graph = CanvaJSLineGraph(csv, csvdata)
            axisX_section = line_graph.set_x_axis({
                'title': 'Title of X data',
                'titleFontSize': 16,
                'labelFontSize': 12,
                'lineColor': '#5BC0DE',
                'titleFontColor': '#5BC0DE',
                'labelFontColor': '#5BC0DE'})
        """
        self._set_axis("axisX", axis_attr)

    def set_axis_y(self, axis_attr=dict()):
        """ Method to configure first axis Y of the line graph.

        :param dict axis_attr: dictionary with canvasjs axisY
            Attributes.

        axisY: http://canvasjs.com/docs/charts/chart-options/axisy/

        Example:

        ::

            line_graph = CanvaJSLineGraph(csv, csvdata)
            line_graph.set_y_axis({'title': 'Title of Y data',
                                   'titleFontSize': 16,
                                   'labelFontSize': 12,
                                   'lineColor': '#5BC0DE',
                                   'titleFontColor': '#5BC0DE',
                                   'labelFontColor': '#5BC0DE'})
        """
        self._set_axis('axisY', axis_attr)

    def set_axis_y2(self, axis_attr=dict()):
        """ Method to configure second axis Y of the line graph.

        :param dict axis_attr: dictionary with canvasjs axisY
            Attributes.

        axisY: http://canvasjs.com/docs/charts/chart-options/axisy/

        Example:

        ::

            line_graph = CanvaJSLineGraph(csv, csvdata)
            axisX_section = line_graph.set_x_axis({
                'title': 'Title of Y data',
                'titleFontSize': 16,
                'labelFontSize': 12,
                'lineColor': '#5BC0DE',
                'titleFontColor': '#5BC0DE',
                'labelFontColor': '#5BC0DE'})
        """
        self._set_axis('axisY2', axis_attr)

    def create_canvasjs(self):
        """ Method to convert all section as javascript function.

        Return a string which contains command line to launch generation of plot,
        js function to create CanvasJS object and the html div that contains
        CanvasJS plot.
        """
        js = self.html_cjs
        # js function to create the chart container
        js_canvas = """
    function drawChart_{0}({1}) {{
        console.log("Start drawChart");
        var chart = new CanvasJS.Chart("chartContainer_{0}",{2}
        );
        chart.render()
    }};
        """
        canvas_attr = super().create_canvas_js_object()
        js += js_canvas.format(self.html_id, self.variables, canvas_attr)
        # js command to run the CanvasJS
        js += """
    $(document).ready(function(){{
        var csv_{0} = document.getElementById('{0}').innerText;
        processData_{0}(csv_{0});
    }});
</script>
        """.format(self.html_id)
        js += self.create_div_chart_container("height: 450px; width: 100%;")
        return js
