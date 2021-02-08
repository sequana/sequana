# -*- coding: utf-8 -*-
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
"""simple summary class to handle summary data with metadata"""
import time
import os
import json
from pathlib import Path

from sequana.lazy import pandas as pd
import colorlog
logger = colorlog.getLogger(__name__)



from sequana.utils.datatables_js import DataTable


__all__ = ["Summary"]



class MultiSummary(object):
    """Helper class to read several json and create summary plots and HTML
content"""

    def __init__(self):
        self.data = {}
        self.order = []

    def read_summary(self, filename, label=None):
        self.filename = filename
        data = json.load(open(self.filename, "r"))
        import os
        if label is None:
            p = Path(filename)
            label = p.name
        self.data[label] = data
        self.order.append(label)

    def remove_summary(self, label):
        if label in self.data and label in self.order:
            del self.data[label] 
            self.order.pop(label)

    def get_html_table(self, user_key_list):
        df = self.get_single_data(user_key_list)
        datatable = DataTable(df, 'name')
        datatable.datatable.datatable_options = {'pageLength': 15,
            'scrollCollapse': 'false',
            'dom': 'Brt',
            'buttons': ['copy', 'csv']}
        js = datatable.create_javascript_function()
        html = datatable.create_datatable(float_format='%.6g')
        return js + html

    def get_single_data(self, user_key_lists):
        # first get the requested data
        data = {}
        for key in self.data.keys():
            values = []
            for user_key in user_key_lists:
                value = self.data[key]
                for depth in user_key.split("/"):
                    value = value[depth]
                values.append(value)
            data[key] = values

        df = pd.DataFrame(data.values(), index=data.keys())
        df.columns = user_key_lists
        df = df.loc[self.order]
        df = df.reset_index() # we need at least one index and one value

        return df 



class Summary(object):
    """

    .. doctest::

        >>> s = Summary("test", "chr1", data={"mean": 1})
        >>> s.name
        sequana_summary_test
        >>> s.sample_name
        chr1


    Here, we prefix the name with the "sequana_summary" tag. Then,
    we populate the sequana version and date automatically. The final
    summary content is then accessible as a dictionary::

        >>> s.as_dict()
        {'data': {'mean': 1},
         'date': 'Thu Jan 18 22:09:13 2018',
         'name': 'sequana_summary_test',
         'sample_name': 'chr1',
         'version': '0.6.3.post1'}

    You can also populate a description dictionary that will provide a
    description for the keys contained in the *data* field. For instance,
    here, the data dictionary contains only one obvious field (mean), we could
    provide a description::

        s.data_description = {"mean": "a dedicated description for the mean"}

    A more general description can also be provided::

        s.description = "bla bla bla"

    """
    def __init__(self, name, sample_name="undefined", data={}, caller=None,
            pipeline_version=None):

        if os.path.exists(name) and name.endswith('json'):
            with open(name, "r") as fin:
                data = json.loads(fin.read())
                self._name = data["name"]
                self.description = data["description"]
                self._data_description = data["data_description"]
                self.sample_name = data["sample_name"]
                self.data = data["data"]
                self.params = data.get("params", {})
                if "caller" in data.keys():
                    self.caller = data["caller"]
                else:
                    self.caller = "undefined"
        else:
            name = name.strip()
            assert len(name.split()) == 1, "no space allowed in the name"
            assert isinstance(data, dict), "data must be a dictionary"
            self._name = name
            self.description = ""
            self._data_description = {}
            self.sample_name = sample_name
            self.data = data
            self.caller = caller
            self.pipeline_version = pipeline_version
            self.params = {}

    def as_dict(self):
        return {
            "name": self.name,
            "sample_name": self.sample_name,
            "version": self.version,
            "pipeline_version": self.pipeline_version,
            "date": self.date,
            "data": self.data,
            "params": self.params,
            "description": self.description,
            "data_description": self.data_description,
            "caller": self.caller,
    }
    def add_params(self, params):
        self.params = params

    def to_json(self, filename):
        import json
        with open(filename, "w") as fh:
            json.dump(self.as_dict(), fh, indent=4, sort_keys=True)

    @property
    def date(self):
        return time.asctime()

    @property
    def name(self):
        return "sequana_summary_" + self._name

    @property
    def version(self):
        from sequana import version
        return version

    @property
    def data_description(self):
        d = {}
        for k in self.data.keys():
            d[k] = self._data_description.get(k, None)
        return d

    @data_description.setter
    def data_description(self, desc):
        self._data_description = {}
        assert isinstance(desc, dict), "data_description must be a dictionary"
        for k,v in desc.items():
            if k not in self.data.keys():
                raise KeyError("{} not a key found in your data dictionary")
            else:
                self._data_description[k] = v
