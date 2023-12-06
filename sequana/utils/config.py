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
""" Sequana report config contains configuration informations to create HTML
report with Jinja2.
"""
import glob
import os
import sys
from datetime import datetime
from pathlib import Path

import easydev
import pkg_resources

time_now = datetime.now().strftime("%m-%d-%Y %H:%M:%S")

# Get sequana informations
version = pkg_resources.get_distribution("sequana").version
script_path = os.path.dirname(os.path.realpath(__file__))

# Default options
output_dir = os.path.realpath(os.getcwd())


def get_entry_name(entry_point):
    return str(entry_point).split("=")[0].strip()


# Modules available
module_dict = {
    get_entry_name(entry_point): entry_point for entry_point in pkg_resources.iter_entry_points("sequana.module")
}

# Find css/js file necessary for the report
sequana_path = easydev.get_package_location("sequana")
css_path = os.sep.join([sequana_path, "sequana", "resources", "css"])
js_path = os.sep.join([sequana_path, "sequana", "resources", "js"])
images_path = Path(sequana_path) / "sequana" / "resources" / "images"

css_list = glob.glob(css_path + os.sep + "*css")
js_list = glob.glob(js_path + os.sep + "*js")
logo = images_path / "logo_256x256.png"

# Sections list for summary.html
summary_sections = list()
