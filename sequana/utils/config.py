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
import inspect
import os
import sys
from datetime import datetime
from pathlib import Path

import sequana
from sequana import version

time_now = datetime.now().strftime("%m-%d-%Y %H:%M:%S")

script_path = os.path.dirname(os.path.realpath(__file__))

# Default options
output_dir = os.path.realpath(os.getcwd())

sequana_path = Path(inspect.getfile(sequana)).parent


css_path = sequana_path / "resources" / "css"
js_path = sequana_path / "resources" / "js"
images_path = sequana_path / "resources" / "images"

css_list = [str(x) for x in css_path.glob("*css")]
js_list = [str(x) for x in js_path.glob("*js")]
logo = images_path / "logo_256x256.png"

# Sections list for summary.html
summary_sections = list()
