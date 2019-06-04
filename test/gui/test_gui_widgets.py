import sys
import os
import time

from sequana.gui import sequana_gui
from sequana import sequana_data, SequanaConfig

import pytest

from tempfile import TemporaryDirectory
from argparse import Namespace




def test_standalone_generic_with_config(qtbot, tmpdir):
    # Standalone for generic case given a wkdir and snakefile

    wkdir = TemporaryDirectory()
    args = Namespace(wkdir=wkdir.name,
                snakefile=sequana_data("test_generic.rules"),
                configfile=sequana_data("test_generic.yaml"))
    widget = sequana_gui.SequanaGUI(ipython=False, user_options=args)
    qtbot.addWidget(widget)
    assert widget.mode == "generic"
    assert widget.generic_factory.is_runnable() == True
    widget.save_project()


    # read back.
    yaml = SequanaConfig(wkdir.name + "/test_generic.yaml").config
    assert yaml['test']["mylist"] == [1,2,3]
