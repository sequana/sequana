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
""" Base class for CanvasJS plot to set legend, title, axis and commun things.
"""


class CanvasJS(object):
    """Base class to create any type of graph with CanvasJS."""

    def __init__(self, html_id):
        """.. rubric:: constructor

        :param str html_id: the ID used in your html. All function
            will have this tag.
        """
        # set id for chart container
        self.html_id = html_id
        # set composite attributes
        self.title_section = dict()
        self.legend_section = dict()
        self.axis_section = dict()
        self.data_section = list()
        # set simple attributes
        self.options = dict()

    def set_title(self, title, title_attr=dict()):
        """Method to configure title of the CanvasJS chart.

        :param str title: title of the chart.
        :param dict title_attr: canvasjs title attributes

        title: http://canvasjs.com/docs/charts/chart-options/title

        Example:

        ::

            canvasjs = CanvasJS()
            canvasjs.set_title("Title of the plot", {'fontFamily': 'verdana',
                                                     'fontSize': 16,
                                                     'verticalAlign': 'top'})

        Populate :attr:`CanvasJS.title_section`
        """
        self.title_section["text"] = title
        self.title_section.update(title_attr)

    def set_legend(self, legend_attr=dict(), hide_on_click=True):
        """Method to configure legend of the CanvasJS chart.

        :param dict legend_attr: dictionary with CanvasJS legend attributes.

        legend: http://canvasjs.com/docs/charts/chart-options/legend/

        Example:

        ::

            canvasjs = CanvasJS()
            canvasjs.set_legend({'fontFamily': 'verdana',
                                 'fontSize': 12,
                                 'verticalAlign': 'bottom'}

        Populate the dictionary :attr:`legend_section.CanvasJS`.
        """
        self.legend_section.update(legend_attr)
        # if you want to hide legend and data on click
        if hide_on_click:
            hide_attr = """function(e){
                if (typeof(e.dataSeries.visible) === "undefined" ||
                    e.dataSeries.visible) {
                    e.dataSeries.visible = false;
                }
                else{
                    e.dataSeries.visible = true;
                }
                chart.render();
            }
            """
            self.legend_section["itemclick"] = hide_attr

    def set_options(self, options=dict()):
        """Method to add options for the CanvasJS chart.

        :param dict options: dictionary with desired options for CanvasJS.

        options: http://canvasjs.com/docs/charts/chart-options/

        Populate the dictionary :attr:`CanvasJS.options`.
        """
        self.options.update(options)

    def _set_axis(self, axis_id, axis_attr=dict()):
        """Method to configure axis.

        :param str axis_id: name of the section. They must have the form of
            this regex "axis[XY][12]".
        :param dict axis_attr: canvasjs axisX/axisY attributes.

        axisX: http://canvasjs.com/docs/charts/chart-options/axisx/
        axisY: http://canvasjs.com/docs/charts/chart-options/axisy/

        Example:

        ::

            canvasjs = CanvasJS()
            axisX_section = canvasjs._set_axis("axisX", {
                'title': 'X num',
                'titleFontSize': 16,
                'labelFontSize': 12,
                'lineColor': '#5BC0DE',
                'titleFontColor': '#5BC0DE',
                'labelFontColor': '#5BC0DE'})

        Populate :attr:`CanvasJS.axis_section`.
        """
        self.axis_section[axis_id] = axis_attr

    def create_canvas_js_object(self):
        """Method to convert all section as javascript function.
        Return string.
        """
        js = "{"
        # create standard options
        if self.options:
            js += self._dict_to_string(self.options) + ","
        # create title and legend sections
        if self.title_section:
            js += "title:{{{0}}},".format(self._dict_to_string(self.title_section))
        if self.legend_section:
            js += "legend:{{{0}}},".format(self._dict_to_string(self.legend_section))
        if self.axis_section:
            axis_list = [
                "{0}:{{{1}}}".format(key, self._dict_to_string(value)) for key, value in self.axis_section.items()
            ]
            js += ",".join(axis_list) + ","
        # dataPoint value must not have quote
        data_list = ["{{{0}}}".format(self._dict_to_string(d)) for d in self.data_section]
        js += "data:[{0}]}}".format(",".join(data_list))
        return js

    def create_div_chart_container(self, style_option=""):
        """HTML div that contains CanvasJS chart.

        :param str style_option: css option for your chart.

        return a string.
        """
        return '<div id="chartContainer_{0}" style="{1}"></div>'.format(self.html_id, style_option)

    def set_data(self, data_dict, index=None):
        """Method to convert dictionnarie as data field for CanvasJS.

        :param dict data_attr: dictionary with CanvasJS data attributes.
        :param int index: which data you want to update.

        data: http://canvasjs.com/docs/charts/chart-options/data/

        Populate :attr:`CanvasJS.data_section`.
        """
        try:
            self.data_section[index].update(data_dict)
        except (IndexError, TypeError):
            self.data_section.append(data_dict)

    def _dict_to_string(self, d):
        """Convert dict to string for CanvasJS.

        Example:

        ::
            dico = {'key1': value1, 'key2': value2, 'key3': value3}
            print(CanvasJS._dict_to_string(dico))

            "key1:value1,key2:value2,key3:value3"
        """
        s = ["{0}:{1}".format(key, self._check_type(value)) for key, value in d.items()]
        return ",\n".join(s)

    def _check_type(self, value):
        """Check type of value to fill javascript sections. String must be
        surrounded by quote and not boolean or integer.

        Javascript variable must not be surrounded by quote. Custom variables
        start with 'data_'.
        """
        try:
            if not value.startswith(("true", "false", "function", "{", "data_")):
                return '"{0}"'.format(value)
        except AttributeError:
            return value
        return value
