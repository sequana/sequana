import sys
import os
import time

from sequana.gui import sequana_gui
from sequana import sequana_data
from sequana import Module

from PyQt5 import QtWidgets as QW
from PyQt5.Qt import QTimer
import pytest

from tempfile import TemporaryDirectory
from argparse import Namespace

from mock import patch
from easydev import TempFile

skiptravis = pytest.mark.skipif("TRAVIS_PYTHON_VERSION" in os.environ, reason="On travis")



@pytest.fixture
def module():
    return Module("quality_control")


def test_settings(qtbot):
    widget = sequana_gui.SequanaGUI(ipython=False)
    qtbot.addWidget(widget)
    widget.read_settings()

    #widget.menuHelp()
    widget.close()


def test_standalone_generic(qtbot, tmpdir, module):

    wkdir = TemporaryDirectory()
    # Standalone for generic case given a wkdir and snakefile (no config)
    args = Namespace(wkdir=wkdir.name,
                snakefile=module.snakefile)
    widget = sequana_gui.SequanaGUI(ipython=False, user_options=args)
    qtbot.addWidget(widget)
    assert widget.mode == "generic"

@pytest.mark.xfail
def test_standalone_generic_with_config(qtbot, tmpdir, module):
    # Standalone for generic case given a wkdir and snakefile
    wkdir = TemporaryDirectory()
    args = Namespace(wkdir=wkdir.name,
                snakefile=module.snakefile, configfile=module.config)
    widget = sequana_gui.SequanaGUI(ipython=False, user_options=args)
    qtbot.addWidget(widget)
    assert widget.mode == "generic"
    assert widget.generic_factory.is_runnable() == True
    widget.save_project()
    widget.unlock_snakemake()
    with TempFile() as fh:
        widget.report_issues(fh.name)

    widget.click_run()

@skiptravis
@pytest.mark.xfail
def test_standalone_generic_with_noconfig_2(qtbot):
    """mimics:

        sequanix -s path_to_snakefile

    followed by selection of a working directory in the GUI + forceall + run
    """
    # From the command line argument
    snakefile = sequana_data("test_gui_generic_snakefile_noconfig.rules")
    wkdir = TemporaryDirectory()
    args = Namespace(snakefile=snakefile)
    widget = sequana_gui.SequanaGUI(ipython=False, user_options=args)
    qtbot.addWidget(widget)
    assert widget.mode == "generic"
    assert widget.generic_factory.is_runnable() is False
    assert widget.ui.run_btn.isEnabled() == False

    widget.generic_factory._directory_browser.set_filenames(wkdir.name)
    widget.force = True
    widget.save_project()

    # check that it worked:
    widget.snakemake_dialog.ui.snakemake_options_general_forceall_value.setChecked(True)
    widget.click_run()
    time.sleep(2)
    widget.click_stop()
    time.sleep(4)
    widget.show_dag()
    time.sleep(4)
    widget.diag.close()

    # in the GUI, we see when it stops. Here, we need to wait a few seconds
    time.sleep(5)
    #data = open(wkdir.name + os.sep + "count.txt").read().split()


def test_open_report(qtbot, tmpdir, module):
    p = tmpdir.mkdir("sub").join('test.html')
    p.write("hello")

    args = Namespace(wkdir=str(p.dirpath()),
                snakefile=module.snakefile)

    # Open a valid HTML file
    widget = sequana_gui.SequanaGUI(ipython=False, user_options=args)
    qtbot.addWidget(widget)
    widget.preferences_dialog.ui.preferences_options_general_htmlpage_value.setText("test.html")
    widget.open_report()
    assert widget.browser.isVisible()
    widget.close()

    # a wrong filename
    widget.preferences_dialog.ui.preferences_options_general_htmlpage_value.setText("dummy.html")
    widget.open_report()
    assert widget.browser.isVisible() is False

    # Now, we unset the working dir and this should just return without openning
    # any report
    widget.sequana_factory._directory_browser.set_empty_path()
    widget.generic_factory._directory_browser.set_empty_path()
    assert widget.working_dir is None
    widget.open_report()
    assert widget.browser.isVisible() is False

    # try using firefox
    widget = sequana_gui.SequanaGUI(ipython=False, user_options=args)
    qtbot.addWidget(widget)
    widget.preferences_dialog.ui.preferences_options_general_htmlpage_value.setText("test.html")
    widget.preferences_dialog.ui.preferences_options_general_browser_value.setCurrentText("firefox")

    def simple_execute(cmd):
        pass
    @patch("easydev.execute", side_effects=simple_execute)
    def runthis(qtbot):
        widget.open_report()

    runthis()
    # There is no dialog browser but firefox should open a tab
    #assert widget.browser.isVisible()


def test_progress_bar(qtbot):
    widget = sequana_gui.SequanaGUI(ipython=False)
    qtbot.addWidget(widget)
    widget.click_run() # defines the regex
    widget.update_progress_bar("0 of 10 steps ")
    widget.update_progress_bar("10 of 10 steps ")
    widget.start_progress()
    widget.end_run()


def test_user_interface_sequana(qtbot):
    widget = sequana_gui.SequanaGUI(ipython=False)
    qtbot.addWidget(widget)
    assert widget.form.count() == 0

    # simulate selection of quality control pipeline
    index = widget.sequana_factory._choice_button.findText("quality_control")
    widget.sequana_factory._choice_button.setCurrentIndex(index)
    widget.ui.tabs_pipeline.setCurrentIndex(0) # set sequana pipeline mode
    widget._update_sequana("quality_control")

    # we should have the focus on the config file now
    assert widget.ui.tabs.currentIndex() == 2
    assert widget.form.count() == 5
    widget.clear_form()
    assert widget.form.count() == 0

    # select no pipeline
    widget._update_sequana('Select a Sequana pipeline')
    assert widget.form.count() == 0


def test_others(qtbot, mocker):
    widget = sequana_gui.SequanaGUI(ipython=True)
    qtbot.addWidget(widget)
    # level and pipeline attribute
    widget.set_level()
    assert widget.sequana_factory.pipeline is None
    # The repr functions
    widget.sequana_factory.__repr__()
    widget.generic_factory.__repr__()

    # menuQuit
    import sequana.gui.messages
    mocker.patch.object(sequana.gui.messages.WarningMessage, "exec_",
        return_value=QW.QMessageBox.Yes)
    widget.menuQuit()

    # menu Help and About
    import sequana.gui.about
    mocker.patch.object(sequana.gui.about.About, "exec_",
        return_value=QW.QMessageBox.Ok)
    widget.menuAbout()
    mocker.patch.object(sequana.gui.help.HelpDialog, "exec_",
        return_value=QW.QMessageBox.Ok)
    widget.menuHelp()


def test_generic_copy_nodir(qtbot):
    # _copy does not work if directory not set
    snakefile = sequana_data("test_gui_generic_snakefile_noconfig.rules")
    configfile = sequana_data("test_gui_generic_config.yaml")
    args = Namespace(snakefile=snakefile, configfile=configfile)
    widget = sequana_gui.SequanaGUI(ipython=False, user_options=args)
    qtbot.addWidget(widget)
    widget.generic_factory._copy_configfile()
    widget.generic_factory._copy_snakefile()


def test_options():
    user_options = sequana_gui.Options()
    options = user_options.parse_args(["--pipeline", "quality_control"])


def test_only(qtbot):
    from easydev import execute
    execute("sequanix --no-splash --testing")


def test_import_config_from_menu(qtbot):
    widget = sequana_gui.SequanaGUI(ipython=True)
    qtbot.addWidget(widget)
    assert widget.sequana_factory._imported_config is None
    # while an existing config file should
    # First, we simulate selection of quality control pipeline
    index = widget.sequana_factory._choice_button.findText("quality_control")
    widget.sequana_factory._choice_button.setCurrentIndex(index)
    widget.ui.tabs_pipeline.setCurrentIndex(0) # set sequana pipeline mode
    widget._update_sequana("quality_control")

    qc = Module("quality_control")
    widget.menuImportConfig(qc.config)
    assert widget.sequana_factory._imported_config is not None

    widget.close()


