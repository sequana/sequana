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
"""Sequana GUI. Can also be used for any snakemake pipeline"""
import sys
import os
import shutil
import re
import time
import psutil
import subprocess as sp
import argparse
import signal
from optparse import OptionParser

from sequana.gui.ui_mainwindow import Ui_MainWindow
from sequana.gui.browser import Browser
from sequana.gui.ipython import QIPythonWidget
from sequana.gui.about import About
from sequana.gui.file_browser import FileBrowser
from sequana.gui.widgets import Ruleform, SVGDialog
from sequana.gui.messages import WarningMessage, CriticalMessage
from sequana.gui.preferences import PreferencesDialog
from sequana.gui.help import HelpDialog
from sequana.gui.snakemake import SnakemakeDialog
from sequana.gui.tools import Tools, QPlainTextEditLogger

from PyQt5 import QtCore, QtGui
from PyQt5 import QtWidgets as QW
from PyQt5.Qt import QTemporaryDir, QMainWindow
from PyQt5.QtCore import Qt

from sequana import snaketools, sequana_data
from sequana.snaketools import Module
from sequana.iotools import YamlDocParser
from sequana import misc

import easydev
import colorlog


def sigint_handler(*args):
    """Handler for the SIGINT signal."""
    sys.stderr.write('\r')
    if QW.QMessageBox.question(None, '', "Are you sure you want to quit?",
                            QW.QMessageBox.Yes | QW.QMessageBox.No,
                            QW.QMessageBox.No) == QW.QMessageBox.Yes:
        QW.QApplication.quit()


class BaseFactory(Tools):
    """Tab on top are all based on this abstract class

    It should provide access to a snakefile and its config file as well
    as working directory.

    Currently, the :class:`SequanaFactory` and :class:`GenericFactory` are
    implemented.

    """
    def __init__(self, mode, run_button):
        super(BaseFactory, self).__init__()
        self.mode = mode
        self._run_button = run_button

        # And finally the working directory
        self._directory_browser = FileBrowser(directory=True)
        self._directory_browser.clicked_connect(self._switch_off_run)

    def _switch_off_run(self):
        self.debug("Switching off run button")
        self._run_button.setEnabled(False)

    def copy(self, source, target, force):
        if os.path.exists(target) and force is False:
            save_msg = WarningMessage(
                    "The file <i>{0}</i> already exists in the working directory".format(source))
            save_msg.setInformativeText(
                    "Do you want to overwrite it?")
            save_msg.setStandardButtons(
                    QW.QMessageBox.Yes | QW.QMessageBox.Discard |
                    QW.QMessageBox.Cancel)
            save_msg.setDefaultButton(QW.QMessageBox.Yes)
            # Yes == 16384
            # Save == 2048
            retval = save_msg.exec_()
            if retval in [16384, 2048]:
                self.warning("Overwritting %s" % target)
                super(BaseFactory, self).copy(source, target)
        else:
                super(BaseFactory, self).copy(source, target)

    def _copy_snakefile(self, force=False):
        if self.snakefile is None:
            self.info("No pipeline selected yet")
            return # nothing to be done

        if self.directory is None:
            self.info("No working directory selected yet (copy snakefile)")
            return

        # source and target filenames
        target = self.directory + os.sep + os.path.basename(self.snakefile)

        if os.path.exists(target) and easydev.md5(target) == easydev.md5(self.snakefile):
            self.info("Target and source (pipeline) are identical. Skipping copy.")
            # if target and source are identical, nothing to do
            return

        # if filename are identical but different, do we want to overwrite it ?
        if os.path.basename(self.snakefile) == target:
            self.warning("%s exists already in %s" % (self.snakefile,
                self.directory))
            return

        self.info("Copying snakefile in %s " % self.directory)
        self.copy(self.snakefile, target, force=force)

    def _copy_configfile(self):
        if self.configfile is None:
            self.info("No config selected yet")
            return # nothing to be done

        if self._directory_browser.path_is_setup() is False:
            self.info("No working directory selected yet (copy config)")
            return

        # FIXME
        # THis does not check the formatting so when saved, it is different
        # from original even though parameters are the same...
        target = self.directory + os.sep + os.path.basename(self.configfile)
        if os.path.exists(target) and easydev.md5(target) == easydev.md5(self.configfile):
            self.info("Target and source (pipeline) are identical. Skipping copy.")
            return

        self.info("Copying config in %s " % self.directory)
        self.copy(self.configfile, self.directory)

    def _get_directory(self):
        filename = self._directory_browser.get_filenames()
        if len(filename):
            return filename
        else:
            return None
    directory = property(_get_directory)

    def __repr__(self):
        return "%s Factory" % self.mode


class SequanaFactory(BaseFactory):
    def __init__(self, run_button, combobox):
        super(SequanaFactory, self).__init__("sequana", run_button)
        self._imported_config = None
        self._choice_button = combobox

        # Some widgets to be used: a file browser for paired files
        fastq_filter = "Fastq file (*.fastq *.fastq.gz *.fq *.fq.gz)"
        self._sequana_paired_tab = FileBrowser(paired=True, file_filter=fastq_filter)
        self._sequana_readtag_label2 = QW.QLabel("Read tag (e.g. _[12].fastq)")
        self._sequana_readtag_lineedit2 = QW.QLineEdit("_R[12]_")

        # Set the file browser input_directory tab
        self._sequana_directory_tab = FileBrowser(directory=True)
        self._sequana_readtag_label = QW.QLabel("Read tag (e.g. _[12].fastq)")
        self._sequana_readtag_lineedit = QW.QLineEdit("_R[12]_")
        self._sequana_pattern_label = QW.QLabel(
            "<div><i>Optional</i> pattern (e.g., Samples_1?/*fastq.gz)</div>")
        self._sequana_pattern_lineedit = QW.QLineEdit()

        # triggers/connectors
        self._sequana_directory_tab.clicked_connect(self._switch_off_run)
        self._choice_button.activated.connect(self._switch_off_run)
        self._sequana_paired_tab.clicked_connect(self._switch_off_run)

    def _get_pipeline(self):
        index = self._choice_button.currentIndex()
        if index == 0:
            return None
        else:
            return self._choice_button.currentText()
    pipeline = property(_get_pipeline)

    def _get_snakefile(self):
        if self.pipeline:
            module = snaketools.Module(self.pipeline)
            return module.snakefile
    snakefile = property(_get_snakefile)

    def _get_configfile(self):
        if self.pipeline:
            module = snaketools.Module(self.pipeline)
            return module.config
    configfile = property(_get_configfile)

    def _get_clusterconfigfile(self):
        if self.pipeline:
            module = snaketools.Module(self.pipeline)
            return module.cluster_config
    clusterconfigfile = property(_get_clusterconfigfile)

    def _get_multiqcconfigfile(self):
        if self.pipeline:
            module = snaketools.Module(self.pipeline)
            return module.multiqc_config
    multiqcconfigfile = property(_get_multiqcconfigfile)

    def _get_schemafile(self):
        if self.pipeline:
            module = snaketools.Module(self.pipeline)
            return module.schema_config
    schemafile = property(_get_schemafile)

    def _get_config(self):
        if self._imported_config:
            cfg = snaketools.SequanaConfig(self._imported_config)
            return cfg
        if self.configfile:
            try:
                cfg = snaketools.SequanaConfig(self.configfile)
                return cfg
            except AssertionError:
                self.warning("Warning: could not parse the config file")
                return
    config = property(_get_config)

    def __repr__(self):
        in1 = self._sequana_directory_tab.get_filenames()
        in2 = self._sequana_paired_tab.get_filenames()
        txt = super(SequanaFactory, self).__repr__()
        txt += "\npipeline:%s\ninput:\n - %s\n - %s\n - directory:%s\n"
        if self.clusterconfigfile:
            txt +=" - cluster config: %s\n" % self.clusterconfigfile
        if self.schemafile:
            txt +=" - schema config: %s" % self.schemafile
        if self.multiqcconfigfile:
            txt +=" - schema config: %s" % self.multiqcconfigfile
        return txt % (self.pipeline, in1, in2, self.directory)


class GenericFactory(BaseFactory):

    def __init__(self, run_button):
        super(GenericFactory, self).__init__("generic", run_button)

        # Define the Snakefile browser and config widgets
        self._snakefile_browser = FileBrowser(directory=False)
        self._config_browser = FileBrowser(directory=False,
            file_filter="YAML file (*.json *.yaml *.yml)")

        # when a snakefile or config is chosen, switch off run button
        self._config_browser.clicked_connect(self._switch_off_run)
        self._snakefile_browser.clicked_connect(self._switch_off_run)
        self._schema = None
        self._multiqcconfigfile = None

    def _return_none(self, this):
        if this is None or len(this) == 0:
            return None
        else:
            return this

    def _get_snakefile(self):
        return self._return_none(self._snakefile_browser.get_filenames())
    snakefile = property(_get_snakefile)

    def _get_configfile(self):
        return self._return_none(self._config_browser.get_filenames())
    configfile = property(_get_configfile)

    def _get_schemafile(self):
        return self._return_none(self._schema)
    schemafile = property(_get_schemafile)

    def _get_multiqcconfigfile(self):
        return self._return_none(self._multiqcconfigfile)
    multiqcconfigfile = property(_get_multiqcconfigfile)

    def _get_config(self):
        filename = self._return_none(self._config_browser.get_filenames())
        if filename:
            try:
                configfile = snaketools.SequanaConfig(filename)
            except AssertionError:
                self.critical("Could not parse the config file %s" % filename)
                return
            except Exception:
                self.critical("Could not parse the config file %s. 2" % filename)
                return
            return configfile
    config = property(_get_config)

    def is_runnable(self):
        flag1 = self._directory_browser.path_is_setup()
        flag2 = self._snakefile_browser.path_is_setup()
        flag3 = self._config_browser.path_is_setup()

        # flag1 and flag2 are compulsary
        # flag3 (configfile) is most tricky to handle since it may be required
        # or not. So we just deal with the flag1 and 2
        if flag1 and flag2:
            return True
        return False

    def __repr__(self):
        txt = super(GenericFactory, self).__repr__()
        txt += "\nsnakefile:%s\nconfigfile:%s\ndirectory:%s\nschema:%s\nmultiqcconfigfile:%s"
        return txt % (self.snakefile, self.configfile, self.directory,
            self.schemafile, self.multiqcconfigfile)


class SequanaGUI(QMainWindow, Tools):
    """

    If quiet, progress bar cannot work.

    - do not copy again requirements if already there
    - extension of the different widgets ?

    Developer Guide
    ------------------

    - The GUI is designed with qt designer as much as possible.
    - All GUI objects are in the **ui** attributes. Additional dialog such as the
      snakemake and preferences dialog have their own modules and stored in attributes
      ending in _dialog

    """
    _not_a_rule = {"requirements", "gatk_bin", "input_directory",
                    "input_pattern", "ignore"}
    _browser_keywords = {"reference"}
    _to_exclude = ["atac-seq", "compressor"]

    def __init__(self, parent=None, ipython=True, user_options={}):
        super(SequanaGUI, self).__init__(parent=parent)

        colorlog.getLogger().setLevel("INFO")
        colorlog.info("Welcome to Sequana GUI (aka Sequanix)")

        self._tempdir = QTemporaryDir()
        self.shell = ""
        self.shell_error = ""
        self._colors = {
            'green': QtGui.QColor(0,170,0),
            'red': QtGui.QColor(170,0,0),
            'orange': QtGui.QColor(170,150,0),
            'blue': QtGui.QColor(0,90,154),
        }

        # some global attributes
        self._undefined_section = "Parameters in no sections/rules"
        #self._config = None

        # Set the regex to catch steps in the progres bar
        self._step_regex = re.compile("([0-9]+) of ([0-9]+) steps")

        self._ipython_tab = ipython
        self.initUI()
        self.read_settings()

        # this should be after initUI and read_settings
        self.set_style_sheet()

        # User option.
        def isset(options, key):
            if key in options and getattr(options,key):
                return True
            else:
                return False

        if isset(user_options, "wkdir"):
            self.info("Setting working directory using user's argument %s" %
                user_options.wkdir)
            if os.path.exists(user_options.wkdir) is False:
                easydev.mkdirs(user_options.wkdir)
            # We must use the absolute path
            abspath = os.path.abspath(user_options.wkdir)
            self.sequana_factory._directory_browser.set_filenames(abspath)
            self.generic_factory._directory_browser.set_filenames(abspath)

        if isset(user_options, "snakefile"):
            filename = user_options.snakefile
            if os.path.exists(filename) is True:
                self.info("Setting snakefile using user's argument %s" %
                    user_options.snakefile)
                self.generic_factory._snakefile_browser.set_filenames(filename)
            else:
                self.error("%s does not exist" % filename)
            self.ui.tabs_pipeline.setCurrentIndex(1)

        if isset(user_options, "configfile"):
            filename = user_options.configfile
            if os.path.exists(filename) is True:
                self.info("Setting config file using user's argument %s" %
                    user_options.configfile)
                self.generic_factory._config_browser.set_filenames(filename)
            self.ui.tabs_pipeline.setCurrentIndex(1)

        if isset(user_options, "pipeline"):
            self.info("Setting Sequana pipeline %s " % user_options.pipeline)
            pipelines = self.sequana_factory.valid_pipelines
            if user_options.pipeline in pipelines:
                index = self.ui.choice_button.findText(user_options.pipeline)
                self.ui.choice_button.setCurrentIndex(index)
                # set focus on pipeline tab
                self.ui.tabs_pipeline.setCurrentIndex(0)
            else:
                self.error("unknown pipeline. Use one of %s " % pipelines)

        if isset(user_options, "input_directory"):
            directory = user_options.input_directory
            self.info("Setting Sequana input directory")
            if directory and os.path.exists(directory) is False:
                self.warning("%s does not exist" % directory)
            elif directory:
                abspath = os.path.abspath(user_options.input_directory)
                self.sequana_factory._sequana_directory_tab.set_filenames(abspath)
            self.ui.tabs_pipeline.setCurrentIndex(0)
            self.ui.tabWidget.setCurrentIndex(0)

        if isset(user_options, "input_files"):
            directory = user_options.input_files
            self.info("Setting Sequana input files")
            dirtab = self.sequana_factory._sequana_paired_tab
            dirtab._set_paired_filenames([os.path.abspath(f)
                for f in user_options.input_files])
            self.ui.tabs_pipeline.setCurrentIndex(0)
            self.ui.tabWidget.setCurrentIndex(1)

        if isset(user_options, "sequana_configfile"):
            cfg = user_options.sequana_configfile
            self.info("Replace Sequana config file")
            self.menuImportConfig(cfg)

        if isset(user_options, "schemafile"):
            schemafile = user_options.schemafile
            self.info("Set the schema file")
            self.menuImportSchema(schemafile)

        # We may have set some pipeline, snakefile, working directory
        self.create_base_form()
        self.fill_until_starting()

    def initUI(self):
        # The logger is not yet set, so we use the module directly
        colorlog.info("Initialising GUI")

        # Set up the user interface from Designer. This is the general layout
        # without dedicated widgets and connections
        self.ui = Ui_MainWindow()
        self.ui.setupUi(self)

        # 2 more dialogs from designer
        self.preferences_dialog = PreferencesDialog(self)
        self.snakemake_dialog = SnakemakeDialog(self)

        self.preferences_dialog.ui.buttonBox.accepted.connect(self.set_level)

        # The IPython dialog, which is very useful for debugging
        if self._ipython_tab is True:
            self.ipyConsole = QIPythonWidget(
                customBanner="Welcome to Sequanix embedded ipython console\n" +
                    "The entire GUI interface is stored in the variable gui\n" +
                    "Note also that you can use this interface as a shell \n" +
                    "command line interface preceding your command with ! character\n")
            #self.ipyConsole.printText("The variable 'foo' andion.")
            self.ipyConsole.execute("from sequana import *")
            self.ipyConsole.execute("import sequana")
            self.ipyConsole.execute("")
            self.ipyConsole.pushVariables({"gui": self})
            self.ui.layout_ipython.addWidget(self.ipyConsole)

        # layout for config file parameters
        widget_form = QW.QWidget()
        self.form = QW.QVBoxLayout(widget_form)
        self.form.setSpacing(10)
        self.ui.scrollArea.setWidget(widget_form)
        self.ui.scrollArea.setWidgetResizable(True)
        self.ui.scrollArea.setMinimumHeight(200)

        # layout for the snakemake output
        self.output = QW.QTextEdit()
        self.output.setReadOnly(True)
        self.ui.layout_snakemake.addWidget(self.output)

        # Add the new logging box widget to the layout
        self.logTextBox = QPlainTextEditLogger(self)
        self.logTextBox.setFormatter(colorlog.ColoredFormatter(
            '%(log_color)s%(asctime)s - %(levelname)s - %(message)s',
            log_colors={
            'DEBUG': 'cyan',
            'INFO': 'green',
            'WARNING': 'yellow',
            'ERROR': 'red',
            'CRITICAL': 'red,bg_white',
            }
            ))
        colorlog.getLogger().addHandler(self.logTextBox)
        self.set_level()
        self.ui.layout_logger.addWidget(self.logTextBox.widget)

        # Connectors to actions related to the menu bar
        self.ui.actionQuit.triggered.connect(self.menuQuit)
        self.ui.action_import_configfile.triggered.connect(self.menuImportConfig)
        self.ui.action_import_schemafile.triggered.connect(self.menuImportSchema)
        self.ui.actionHelp.triggered.connect(self.menuHelp)
        self.ui.actionAbout.triggered.connect(self.menuAbout)
        self.ui.actionSnakemake.triggered.connect(self.snakemake_dialog.exec_)
        self.ui.actionPreferences.triggered.connect(self.preferences_dialog.exec_)

        self.preferences_dialog.ui.preferences_options_general_tooltip_value.clicked.connect(self.set_style_sheet)

        # connectors related to the pipeline tabs (pipeline/generic)
        self.set_sequana_pipeline()
        self.set_generic_pipeline()

        # The run/save/dag footer buttons
        self.connect_footer_buttons()

        self.process = QtCore.QProcess(self)
        self.process.started.connect(lambda: self.ui.run_btn.setEnabled(False))
        self.process.started.connect(lambda: self.ui.stop_btn.setEnabled(True))
        self.process.started.connect(lambda: self.ui.unlock_btn.setEnabled(False))
        self.process.started.connect(lambda: self.start_progress)
        self.process.started.connect(lambda: self.ui.save_btn.setEnabled(False))
        self.process.started.connect(lambda: self.ui.tabs_pipeline.setEnabled(False))

        self.process.finished.connect(lambda: self.ui.run_btn.setEnabled(True))
        self.process.finished.connect(lambda: self.ui.stop_btn.setEnabled(False))
        self.process.finished.connect(lambda: self.ui.unlock_btn.setEnabled(True))
        self.process.finished.connect(lambda: self.ui.save_btn.setEnabled(True))
        self.process.finished.connect(lambda:self.ui.tabs_pipeline.setEnabled(True))
        self.process.finished.connect(self.end_run)

        self.process.readyReadStandardOutput.connect(
            self.snakemake_data_stdout)
        self.process.readyReadStandardError.connect(self.snakemake_data_error)

        # This is for the show dag btn. Created here once for all
        self.process1 = QtCore.QProcess(self)
        self.process2 = QtCore.QProcess(self)

        self.ui.tabWidget.currentChanged.connect(lambda: self.ui.run_btn.setEnabled(False))

        # if we are on one of those clusters, switch to the cluster choice in
        # the pipeline control combo box
        if misc.on_cluster(['tars-']) is True:
             self.ui.comboBox_local.setCurrentText("cluster")

        # connect show advanced button with the until/starting frame
        self.ui.show_advanced_control.clicked.connect(self.click_advanced)
        self.ui.frame_control.hide()

    def _get_opacity(self):
        dialog = self.preferences_dialog
        box = dialog.ui.preferences_options_general_tooltip_value
        if box.isChecked():
            return 255
        else:
            return 0
    tooltip_opacity = property(_get_opacity)

    def set_style_sheet(self):
        self.setStyleSheet("""QToolTip {
                           background-color: #aabbcc;
                           color: black;
                           border-style: double;
                           border-width: 3px;
                           border-color: green;
                           border-radius: 5px;
                           margin:5px;
                           opacity: %s;
                           } ;

                            """ % self.tooltip_opacity)


    #|-----------------------------------------------------|
    #|                       MENU related                  |
    #|-----------------------------------------------------|
    def menuImportConfig(self, configfile=None):
        # The connector send a False signal but default is None
        # so we need to handle the two cases
        if self.snakefile is None:
            self.error("You must set a pipeline first")
            msg = WarningMessage(("You must set a pipeline first"))
            msg.exec_()
            return

        if configfile and os.path.exists(configfile) is False:
            self.error("Config file (%s) does not exists" % configfile)
            return

        if configfile is None or configfile is False:
            self.info("Importing config file.")
            file_filter="YAML file (*.json *.yaml *.yml)"
            browser = FileBrowser(file_filter=file_filter)
            browser.browse_file()
            configfile = browser.paths

        if configfile:
            self.sequana_factory._imported_config = configfile
        else:
            self.sequana_factory._imported_config = None
        self.create_base_form()

    def menuImportSchema(self, schemafile=None ):
        if schemafile:
            self.generic_factory._schema = schemafile
            return

        self.info("Importing YAML schema file.")
        file_filter="YAML file (*.yaml *.yml)"
        browser = FileBrowser(file_filter=file_filter)
        browser.browse_file()
        schemafile = browser.paths
        if schemafile:
            self.generic_factory._schema = schemafile
        else:
            self.generic_factory._schema = None

    def menuAbout(self):
        from sequana import version
        url = 'sequana.readthedocs.io'
        widget = About()
        widget.setText("Sequana version %s " % version)
        widget.setInformativeText("""
            Online documentation on <a href="http://%(url)s">%(url)s</a>
            <br>
            <br>
            Authors: Thomas Cokelaer and Dimitri Desvillechabrol, 2017-2018
            """ % {"url": url})
        widget.setWindowTitle("Sequana")
        #widget.setStandardButtons(QW.QMessageBox.Ok)
        retval = widget.exec_()
        if retval == QW.QMessageBox.Ok:
            widget.close()

    def menuHelp(self):
        url = 'sequana.readthedocs.io'
        pipelines_text = "<ul>\n"
        url = "http://sequana.readthedocs.io/en/master"
        for pipeline in snaketools.pipeline_names:
            pipelines_text += '    <li><a href="%(url)s/pipeline_%(name)s.html">%(name)s</a></li>\n' %\
              {"url":url,"name":pipeline}
        pipelines_text += "</ul>"

        msg = HelpDialog(pipelines=pipelines_text)
        retval = msg.exec_()
        if retval == QW.QMessageBox.Ok:
            msg.close()

    def menuQuit(self):
        self._quit_msg = WarningMessage("Do you really want to quit ?")
        self._quit_msg.setStandardButtons(QW.QMessageBox.Yes | QW.QMessageBox.No)
        self._quit_msg.setDefaultButton(QW.QMessageBox.No)
        quit_answer = self._quit_msg.exec_()
        if quit_answer == QW.QMessageBox.Yes:
            self.close()

    def set_level(self):
        # Set the level of the logging system
        pref = self.preferences_dialog.ui
        level = pref.preferences_options_general_logging_value.currentText()
        level =  getattr(colorlog.logging.logging, level)
        colorlog.getLogger().setLevel(level)

    # ---------------------------------------------------------------
    # More GUI / reading the snakefile (sequana or generic)
    # ---------------------------------------------------------------
    def set_sequana_pipeline(self):
        # The pipeline connectors
        pipelines = sorted(snaketools.pipeline_names)
        pipelines = [this for this in pipelines if this not in self._to_exclude]
        self.ui.choice_button.addItems(pipelines)
        self.ui.choice_button.activated[str].connect(self._update_sequana)

        # FIXME do we want to use this ?
        self.ui.choice_button.installEventFilter(self)

        # populate the factory with the choice button
        self.sequana_factory = SequanaFactory(
            combobox=self.ui.choice_button,
            run_button=self.ui.run_btn)
        self.sequana_factory.valid_pipelines = pipelines

        # a local alias
        saf = self.sequana_factory

        # add widgets for the working dir
        self.ui.layout_sequana_wkdir.addWidget(saf._directory_browser)

        # add widget for the input sample
        #self.ui.layout_sequana_input_files.addWidget(saf._sequana_paired_tab)
        #hlayout = QW.QHBoxLayout()
        #hlayout.addWidget(saf._sequana_readtag_label2)
        #hlayout.addWidget(saf._sequana_readtag_lineedit2)
        #self.ui.layout_sequana_input_files.addLayout(hlayout)

        # add widget for the input directory
        self.ui.layout_sequana_input_dir.addWidget(saf._sequana_directory_tab)
        hlayout = QW.QHBoxLayout()
        hlayout.addWidget(saf._sequana_readtag_label)
        hlayout.addWidget(saf._sequana_readtag_lineedit)
        self.ui.layout_sequana_input_dir.addLayout(hlayout)
        hlayout = QW.QHBoxLayout()
        hlayout.addWidget(saf._sequana_pattern_label)
        hlayout.addWidget(saf._sequana_pattern_lineedit)
        self.ui.layout_sequana_input_dir.addLayout(hlayout)


    @QtCore.pyqtSlot(str)
    def _update_sequana(self, index):
        """ Change options form when user changes the pipeline."""
        if self.ui.choice_button.findText(index) == 0:
            self.clear_form()
            self.rule_list = []
            self.fill_until_starting()
            return

        self.info("Reading sequana %s pipeline" % index)
        self.create_base_form()
        # Is there a cluster config file ?
        dialog = self.snakemake_dialog.ui

        if self.sequana_factory.clusterconfigfile:
            dialog.snakemake_options_cluster_cluster__config_value.set_filenames(
                self.sequana_factory.clusterconfigfile)
        else:
            dialog.snakemake_options_cluster_cluster__config_value.set_filenames("")
        self.fill_until_starting()
        self.switch_off()
        # Reset imported config file in SequanaFactory
        self.sequana_factory._imported_config = None

    def set_generic_pipeline(self):

        self.generic_factory = GenericFactory(self.ui.run_btn)
        gaf = self.generic_factory

        # The config file connectors
        gaf._config_browser.clicked_connect(self.create_base_form)

        # Update the main UI with
        self.ui.layout_generic_snakefile.addWidget(gaf._snakefile_browser)
        self.ui.layout_generic_config.addWidget(gaf._config_browser)
        self.ui.layout_generic_wkdir.addWidget(gaf._directory_browser)

        # When user press the cancel button, the config file browser is reset
        self.ui.cancel_push_button.clicked.connect(
            self.generic_factory._config_browser.set_empty_path)

    # ---------------------------------------------------------------------
    # Footer connectors
    # ---------------------------------------------------------------------

    def connect_footer_buttons(self):
        self.ui.run_btn.setEnabled(False)
        self.ui.run_btn.clicked.connect(self.click_run)

        self.ui.stop_btn.clicked.connect(self.click_stop)
        self.ui.stop_btn.setEnabled(False)

        self.ui.unlock_btn.clicked.connect(self.ui.run_btn.setEnabled)
        self.ui.unlock_btn.clicked.connect(self.unlock_snakemake)
        self.ui.unlock_btn.setEnabled(True)

        self.ui.report_btn.setEnabled(True)
        self.ui.report_btn.clicked.connect(self.open_report)

        self.ui.save_btn.clicked.connect(self.save_project)

        self.ui.dag_btn.setEnabled(False)
        self.ui.dag_btn.clicked.connect(self.show_dag)

    # -----------------------------------------------------------------
    # function to link to the factory (sequana or generic)
    # -----------------------------------------------------------------

    def _get_mode(self):
        # figure out if we are dealing with a sequana pipeline or not
        index = self.ui.tabs_pipeline.currentIndex()
        if index == 0:
            return "sequana"
        elif index == 1:
            return "generic"
    mode = property(_get_mode)

    def _get_factory(self):
        return getattr(self, "%s_factory" % self.mode)
    factory = property(_get_factory)

    def _get_config(self):
        return getattr(self, "%s_factory" % self.mode).config
    config = property(_get_config)

    def _get_configfile(self):
        return getattr(self, "%s_factory" % self.mode).configfile
    configfile = property(_get_configfile)

    def _get_snakefile(self):
        return getattr(self, "%s_factory" % self.mode).snakefile
    snakefile = property(_get_snakefile)

    def _get_working_dir(self):
        return getattr(self, "%s_factory" % self.mode).directory
    working_dir = property(_get_working_dir)

    # ----------------------------------------------------------------------
    # Snakemake related (config, running)
    # ----------------------------------------------------------------------

    def fill_until_starting(self):
        active_list = [w.get_name() for w in self.rule_list if w.get_do_rule()]

        self.ui.until_box.clear()
        self.ui.until_box.addItems([None] + active_list)

        self.ui.starting_box.clear()
        self.ui.starting_box.addItems([None] + active_list)

    # ----------------------------------------------------------
    #  Config file related
    # ---------------------------------------------------------

    def _set_focus_on_config_tab(self):
        # Set focus on config file
        if self._ipython_tab:
            self.ui.tabs.setCurrentIndex(3)
        else:
            self.ui.tabs.setCurrentIndex(2)

    # --------------------------------------------------------------------
    #  Advanced control
    # --------------------------------------------------------------------
    def click_advanced(self):
        if self.ui.frame_control.isHidden():
            self.ui.frame_control.show()
        else:
            self.ui.frame_control.hide()

    # --------------------------------------------------------------------
    # Others
    # --------------------------------------------------------------------

    def clear_layout(self, layout):
        """ Clean all widgets contained in a layout. """
        while layout.count():
            child = layout.takeAt(0)
            if child.widget() is not None:
                child.widget().deleteLater()
            elif child.layout() is not None:
                self.clear_layout(child.layout())

    # --------------------------------------------------------------------
    # Running snakemake
    # --------------------------------------------------------------------

    def _clean_line(self, line):
        # TODO: surely there is a better way to do that and not overlap
        # with tools.py ...
        line = line.replace("b'\\r'", "")
        line = line.replace("b'\r'", "")
        line = line.replace("b'\\r '", "")
        line = line.replace("b'\r '", "")
        line = line.replace("b' '", "")
        line = line.replace("\\t", "&nbsp;"*4)
        line = line.replace("'b'", "")
        for this in ["b'", 'b"', "\r"]:
            if line.startswith(this):
                line = line.replace(this, "")
        if line.startswith('b"'):
            line = line.replace('b"', "")
        line = line.rstrip("\\x1b[0m")
        line = line.replace("\\x1b[33m", "")
        return line

    def snakemake_data_stdout(self):
        """ Read standard output of snakemake process """
        data = str(self.process.readAllStandardOutput())
        self.shell += data
        self.update_progress_bar(data)

        for this in data.split("\\n"):
            line = this.strip()
            if line and len(line) > 3 and "complete in" not in line: # prevent all b'' strings
                line = self._clean_line(line)
                if len(line.strip()) == 0:
                    continue
                self.output.append('<font style="color:blue">' + line +'</font>')

    def snakemake_data_error(self):
        """ Read error output of snakemake process """
        error = str(self.process.readAllStandardError())
        self.shell_error += error
        self.update_progress_bar(error)
        for this in error.split("\\n"):
            line = this.strip()
            if line and len(line) > 3 and "complete in" not in line: # prevent all b'' strings
                line = self._clean_line(line)
                if line.startswith("b'"):
                    line = line[2:]
                    line.rstrip("'")
                grouprex = self._step_regex.findall(line)
                if grouprex:
                    self.output.append('<font style="color:orange">' + line +'</font>')
                elif "Error" in line:
                    self.output.append('<font style="color:red">' + line +'</font>')
                else:
                    self.output.append('<font style="color:green">' + line +'</font>')

    def get_until_starting_option(self):
        """ Return list with starting rule and end rule.
        """
        until_rule = self.ui.until_box.currentText()
        starting_rule = self.ui.starting_box.currentText()
        option = []
        if until_rule:
            option += ["--no-hooks", "-U", until_rule]
        if starting_rule:
            option += ["-R", starting_rule]
        return option

    def _get_snakemake_command(self, snakefile):
        """If the cluster option is selected, then the cluster field in
        the snakemake menu must be set to a string different from empty string.

        If we are on TARS, we also must set the option to cluster (not local)

        If one of the previous cases is true, this function returns None

        """
        dialog = self.snakemake_dialog      # an alias
        snakemake_line = ["-s", snakefile, "--stat", "stats.txt", "-p"]

        if self.ui.comboBox_local.currentText() == "local":
            if misc.on_cluster(["tars-"]):
                msg = WarningMessage(("You are on TARS cluster. Please set the"
                    "batch options and select the cluster option (not local)"))
                msg.exec_()
                return None
            snakemake_line += dialog.get_snakemake_local_options()
        elif self.ui.comboBox_local.currentText() == "cluster":
            cluster = dialog.ui.snakemake_options_cluster_cluster_value.text()
            if len(cluster.strip()) == 0:
                msg = WarningMessage(("You selected a 'cluster run' but the "
                    "cluster preferences are not set. Either switch to a local "
                    "run or set a correct string in the Snakemake options menu "
                    "(in cluster tab/ cluster field.)"))
                msg.exec_()
                return None
            snakemake_line += dialog.get_snakemake_cluster_options()

            #cluster_config = dialog.ui.snakemake_options_cluster_config_value.text()
            #cluster_config = cluster_config.strip()
            #if len(cluster_config):
            #    snakemake_line += ["--cluster-config", cluster_config]

        snakemake_line += dialog.get_snakemake_general_options()
        snakemake_line += self.get_until_starting_option()

        others = self.snakemake_dialog.ui.snakemake_options_general_custom.text()
        if others.strip():
            snakemake_line += others.split()

        if self.configfile:
            configfile = os.path.basename(self.configfile)
            snakemake_line += ["--configfile", configfile]

        return snakemake_line

    def _set_pb_color(self, color):
        self.ui.progressBar.setStyleSheet("""
            QProgressBar {{
                color: black;
                border: 2px solid grey;
                margin: 2px;
                border-radius: 5px;
                text-align: center;
            }}
            QProgressBar::chunk {{
                background: {};
                }}""".format(color))
        #pal = self.ui.progressBar.palette()
        #pal.setColor(QtGui.QPalette.Highlight, self._colors['blue'])
        #self.ui.progressBar.setPalette(pal)

    def click_run(self):
        # set focus on the snakemake output
        if self.snakefile is None or self.working_dir is None:
            self.warning('Working directory or snakefile not set.')
            return
        self.ui.tabs.setCurrentIndex(0)
        self.shell_error = ""
        self.shell = ""

        # Prepare the command and working directory.
        if self.working_dir is None:
            self.warning("Set the working directory first")
            return

        # We copy the sequana and genereic snakefile into a filename called
        # Snakefile
        snakefile = self.working_dir + os.sep + os.path.basename(self.snakefile)

        if os.path.exists(snakefile) is False:
            self.critical("%s does not exist" % snakefile)
            return

        snakemake_args = self._get_snakemake_command(snakefile)
        if snakemake_args is None:
            return

        # the progress bar
        self._set_pb_color(self._colors['blue'].name())
        self.ui.progressBar.setValue(1)

        # Start process
        # If an argument contains spaces, we should use quotes. However,
        # with PyQt quotes must be escaped

        args = []
        for this in snakemake_args:
            if re.search(r"\s", this) is True:
                args.append("\"%s\"" % this)
            else:
                args.append(this)
        snakemake_args = args
        self.info("Starting process with snakemake %s " % " ".join(snakemake_args))
        self.output.clear()
        self.process.setWorkingDirectory(self.working_dir)
        self.process.start("snakemake", snakemake_args)

    # -------------------------------------------------------------------
    # Create the base form
    # -------------------------------------------------------------------
    def create_base_form(self):
        """ Create form with all options necessary for a pipeline.

        ::

            ########################################################
            #   valid python docstring to be interepreted by sphinx
            #
            #   section:
            #      item1: 10
            #      item2: 20

        """
        self.rule_list = []
        if self.config is None:
            self.clear_form()
            return

        self.info("Creating form based on config file")
        self.clear_form()
        rules_list = list(self.config._yaml_code.keys())

        # We do not sort the list of rules anymore so that it is like in the
        # config file
        #rules_list.sort()
        self.necessary_dict = {}

        # For each section, we create a widget (RuleForm). For isntance, first,
        # one is accessible as follows: gui.form.itemAt(0).widget()

        docparser = YamlDocParser(self.configfile)
        import ruamel.yaml.comments
        for count, rule in enumerate(rules_list):
            self.debug("Scanning rule %s" % rule)
            # Check if this is a dictionnary
            contains = self.config._yaml_code[rule]

            # If this is a section/dictionary, create a section

            if isinstance(contains, (ruamel.yaml.comments.CommentedMap, dict)) and (
                    rule not in SequanaGUI._not_a_rule):
                # Get the docstring from the Yaml section/rule
                docstring = docparser._block2docstring(rule)

                # Get any special keywords
                specials = docparser._get_specials(rule)

                #self.ui.preferences_options_general_addbrowser_value
                dialog = self.preferences_dialog.ui
                option = dialog.preferences_options_general_addbrowser_value.text()
                option = option.strip()
                option = option.replace(";", " ").replace(",", " ")

                if len(option):
                    keywords = option.split()
                else:
                    keywords = []
                keywords += self._browser_keywords
                keywords = list(set(keywords))

                rule_box = Ruleform(rule, contains, count, keywords, specials=specials)
                rule_box.connect_all_option(
                    lambda: self.ui.run_btn.setEnabled(False))

                # Try to interpret it with sphinx
                from sequana.misc import rest2html
                try:
                    self.debug("parsing docstring of %s" % rule)
                    comments = rest2html(docstring).decode()
                    rule_box.setToolTip(comments)
                except Exception as err:
                    print(err)
                    self.warning("Could not interpret docstring of %s" % rule)
                    rule_box.setToolTip("")

                self.form.addWidget(rule_box)
                self.rule_list.append(rule_box)
                rule_box.connect_do(self.fill_until_starting)
            else:
                # this is a parameter in a section, which may be
                # a list, a None or something else
                if isinstance(contains, list):
                    self.necessary_dict = dict(self.necessary_dict,
                                           **{rule: contains})
                elif contains is None or contains in ["''", '""']:
                    self.necessary_dict = dict(self.necessary_dict,
                                           **{rule: None})
                else:
                    self.necessary_dict = dict(self.necessary_dict,
                                           **{rule: '{0}'.format(contains)})

        # if this is a generic pipeline, you may have parameters outside of a
        # section
        if self.mode == "generic" and len(self.necessary_dict):
            rule_box = Ruleform(self._undefined_section,
                                self.necessary_dict, -1, generic=True)
            self.form.addWidget(rule_box)
        self._set_focus_on_config_tab()

    # ----------------------------------------------------------
    # STOP footer button
    # ----------------------------------------------------------

    def click_stop(self):
        """The stop button"""
        self._set_pb_color(self._colors['orange'].name())

        # For windows:
        # http://stackoverflow.com/questions/8232544/how-to-terminate-a-process-without-os-kill-osgeo4w-python-2-5

        if self.process.state() != 0:
            pid = self.process.pid()
            self.warning("Process {} running , stopping it... ".format(pid))
            # We must use a ctrl+C interruption so that snakemake
            # handles the interruption smoothly. However, child processes
            # are lost so we also need to get their IDs and kill them.
            self.info("killing the main snakemake process. This may take a few seconds ")
            try:
                self.info("process pid={} being killed".format(self.process.pid()))
                pid_children = [this.pid for this in
                                psutil.Process(pid).children(recursive=True)]
                # Kills the main process
                os.kill(pid, signal.SIGINT)
                # And the children
                for this in pid_children:
                    self.info("Remove pid {} ".format(this))
                    try:
                        os.kill(this, signal.SIGINT)
                    except Exception as err:
                        print(err)
                time.sleep(4)

            except Exception as err:
                print(err)
                pass # already stopped ?
            self.info("Process killed successfully.")
        self.ui.save_btn.setEnabled(True)
        self.ui.run_btn.setEnabled(True)
        self.ui.stop_btn.setEnabled(False)
        self.ui.tabs_pipeline.setEnabled(True)

    # --------------------------------------------------------------------
    # Progress bar
    # --------------------------------------------------------------------

    def update_progress_bar(self, line):
        """ Parse with a regex to retrieve current step and total step.
        """
        grouprex = self._step_regex.findall(line)
        # Use last "x of y" (not the first item at position 0)
        if grouprex:
            step = int(grouprex[-1][0]) / float(grouprex[-1][1]) * 100
            self.ui.progressBar.setValue(step)
        if "Nothing to be done" in line:
            self.ui.progressBar.setValue(100)

    def start_progress(self):
        self.ui.progressBar.setRange(0,1)

    def end_run(self):
        pal = self.ui.progressBar.palette()
        if self.ui.progressBar.value() >= 100 :
            self._set_pb_color(self._colors['green'].name())
            self.info('Run done. Status: successful')
        else:
            self._set_pb_color(self._colors['red'].name())
            text = 'Run manually to check the exact error or check the log.'
            if "--unlock" in self.shell_error:
                text += "<br>You may need to unlock the directory. "
                text += "click on Unlock button"
                self.critical(text)
            return

    def _get_force(self):
        dialog = self.preferences_dialog
        box = dialog.ui.preferences_options_general_overwrite_value
        return box.isChecked()
    def _set_force(self, boolean):
        assert boolean in [True, False]
        dialog = self.preferences_dialog
        box = dialog.ui.preferences_options_general_overwrite_value
        box.setChecked(boolean)
    force = property(_get_force, _set_force)

    def save_project(self):
        self.info('Saving project')

        if self.configfile is None:
            if self.mode == "generic":
                if self.generic_factory.is_runnable():
                    self.critical("save_project: Generic case without config file")
                    self._save_teardown()
                else:
                    msg = WarningMessage("You must select a Snakefile and a working directory.")
                    msg.exec_()
            elif self.mode == "sequana":
                msg = WarningMessage("You must choose a pipeline first.")
                msg.exec_()
            return

        if self.working_dir is None:
            self.critical('save_project: no working dir: return')
            msg = WarningMessage("You must select a working directory first.")
            msg.exec_()
            return

        try:
            form_dict = dict(self.create_form_dict(self.form),
                                 **self.necessary_dict)
        except AttributeError as err:
            self.error(err)
            msg = WarningMessage("You must choose a pipeline before saving.")
            msg.exec_()
            return

        # get samples names or input_directory
        if self.mode == "sequana":
            self.info("Sequana case")
            flag1 = self.sequana_factory._sequana_directory_tab.get_filenames()
            flag2 = self.sequana_factory._sequana_paired_tab.get_filenames()

            if self.ui.tabWidget.currentIndex() == 0 and len(flag1) == 0 or \
                    self.ui.tabWidget.currentIndex() == 1 and len(flag2) == 0:
                msg = WarningMessage("You must choose an input first.")
                msg.exec_()
                return

            filename = self.sequana_factory._sequana_directory_tab.get_filenames()
            form_dict["input_directory"] = (filename)

            # If pattern provided, the input_directory is reset but used in
            # the pattern as the basename
            pattern = self.sequana_factory._sequana_pattern_lineedit.text()
            if len(pattern.strip()):
                form_dict["input_pattern"] = (filename)
                form_dict["input_pattern"] += os.sep + pattern.strip()
                form_dict["input_directory"] = ""

            readtag = self.sequana_factory._sequana_readtag_lineedit.text()
            if len(readtag.strip()):
                form_dict["input_readtag"] = readtag
            else:
                form_dict["input_readtag"] = "_R[12]_"

        elif self.mode == "generic":
            # Here we save the undefined section in the form.
            if self._undefined_section in form_dict.keys():
                for key, value in form_dict[self._undefined_section].items():
                    form_dict[key] = value
                del form_dict[self._undefined_section]
                self.info("Generic case")

        # Let us update the attribute with the content of the form
        # This uses the user's information

        cfg = self.config
        cfg.config.update(form_dict)
        cfg._update_yaml()
        cfg.cleanup()
        self.cfg = cfg

        pref = self.preferences_dialog.ui
        box = pref.preferences_options_general_schema_value
        checked_schema = box.isChecked()

        if self.working_dir:
            # Save the configuration file
            if self.mode == "sequana":
                yaml_path = self.working_dir + os.sep + "config.yaml"
                self.warning("copy requirements (if any)")
                cfg.copy_requirements(target=self.working_dir)
            elif self.mode == "generic":
                yaml_path = self.working_dir + os.sep + \
                                os.path.basename(self.generic_factory.configfile)

            if os.path.isfile(yaml_path) and self.force is False:
                save_msg = WarningMessage(
                        "The file <i>{0}</i> already exist".format(yaml_path))
                save_msg.setInformativeText(
                        "Do you want to overwrite the file?")
                save_msg.setStandardButtons(
                        QW.QMessageBox.Yes | QW.QMessageBox.Discard |
                        QW.QMessageBox.Cancel)
                save_msg.setDefaultButton(QW.QMessageBox.Yes)
                # Yes == 16384
                # Save == 2048
                retval = save_msg.exec_()
                if retval in [16384, 2048]:
                    self.info("Saving config file (exist already)")
                    if checked_schema is False:
                        cfg.save(yaml_path, cleanup=False)
                    else:
                        ret = self._check_and_save_config(cfg, yaml_path)
                        if ret is False:
                            # we do not want to save the config file and call
                            # _save_teardown
                            return
            else:
                self.warning("Saving config file (does not exist)")
                if checked_schema is False:
                    cfg.save(yaml_path, cleanup=False)
                else:
                    ret = self._check_and_save_config(cfg, yaml_path)
                    if ret is False:
                        # we do not want to save the config file and call
                        # _save_teardown
                        return

            # Save the configuration file for the cluster
            if self.mode == "sequana" and self.sequana_factory.clusterconfigfile:
                target = self.working_dir + os.sep + "cluster_config.json"
                shutil.copy(self.sequana_factory.clusterconfigfile, target)
                # replace the name of the original file with the target one so
                # that the target can be edited. The target will also be used in
                # place of the original version when launnching snakemake!
                #self.snakemake_dialog.set()
                self.snakemake_dialog.ui.snakemake_options_cluster_cluster__config_value.set_filenames(target)

            # Save the multiqc_config file if provided in sequana pipeline
            if self.mode == "sequana" and self.sequana_factory.multiqcconfigfile:
                target = self.working_dir + os.sep + "multiqc_config.yaml"
                shutil.copy(self.sequana_factory.multiqcconfigfile, target)


        else:
            self.critical("Config file not saved (no wkdir)")
            msg = WarningMessage("You must set a working directory", self)
            msg.exec_()
            self.switch_off()
            return

        self._save_teardown()

    def _save_teardown(self):
        # Finally, save project and update footer run button
        self.factory._copy_snakefile(self.force)
        self.debug('Switching RUN and DAG button on')
        self.ui.run_btn.setEnabled(True)
        self.ui.dag_btn.setEnabled(True)

    def _check_and_save_config(self, cfg, yaml_path):
        # Here we save the config.yaml file when changed
        # However, before that if there is a schema, we can 
        # use it. This is the case for some sequana pipelines

        # return False if the config is invalid and do not save it

        if self.mode == "sequana" and \
                self.sequana_factory.schemafile is None:
            self.warning("No Schema found to validate the config file")

        if self.mode == "sequana" and self.sequana_factory.schemafile:
            schemafile = self.sequana_factory.schemafile
        elif self.mode == "generic" and self.generic_factory.schemafile:
            schemafile = self.generic_factory.schemafile
        else:
            schemafile = None

        if schemafile:
            # check that the config file is correct before saving it
            # only if we have a schema_config file.
            self.info("Checking config file with provided schema file.")
            # We save the config as a dummy temporary file to check it
            # if correct, we then save the file. If not, we provide an help
            # message
            from easydev import TempFile
            with TempFile(suffix=".yaml") as fout:
                # save a temporary version
                cfg.save(fout.name, cleanup=False)
                import ruamel
                import warnings
                from pykwalify.core import Core
                # causes issue with ruamel.yaml 0.12.13. Works for 0.15
                try:warnings.simplefilter('ignore', ruamel.yaml.error.UnsafeLoaderWarning)
                except: pass

                try:
                    # open the config and the schema file
                    c = Core(source_file=fout.name, schema_files=
                        [schemafile])
                except Exception as err:
                    print(err)
                    return False

                try:
                    c.validate()
                except Exception as err:
                    print(err)
                    error_msg = "<b>CRITICAL: INVALID CONFIGURATION FILE</b>\n"
                    error_msg += "<pre>" + str(err) + "</pre>"
                    self.critical(error_msg)
                    self.switch_off()
                    msg = WarningMessage(error_msg, self)
                    msg.exec_()
                    return False
        cfg.save(yaml_path, cleanup=False)

    def switch_off(self):
        self.debug('Switching RUN and DAG button off')
        self.ui.run_btn.setEnabled(False)
        self.ui.dag_btn.setEnabled(False)

    def _reset_schema(self):
        self.schemafile = None
    # -----------------------------------------------------------------------
    # SAVE LOG in a files
    # -----------------------------------------------------------------------

    def report_issues(self, filename="issue_debug.txt"):
        # save shell + shell_error in working directory as well as snakemake and
        # config file.
        with open(filename, "w") as fh:
            fh.write("\nsequanix logger  ----------------------------------\n")
            try:
                file_logger = self.save_logger()
                with open(file_logger, "r") as fin:
                    fh.write(fin.read())
            except: pass

            fh.write("\nsequanix shell   ----------------------------------\n")
            try:fh.writelines(self.shell)
            except:fh.write("No shell info")

            fh.write("\nsequanix shell error ------------------------------\n")
            try:fh.writelines(self.shell_error)
            except:fh.write("No shell error info")
        url = "https://github.com/sequana/sequana/issues "
        print("Created a file called {} to be posted on {}.".format(
            filename, url))
        self.init_logger()

    # -----------------------------------------------------------------------
    # UNLOCK footer button
    # -----------------------------------------------------------------------

    def unlock_snakemake(self):
        if self.working_dir is None or self.snakefile is None:
            self.warning("working directory or snakefile not set")
            return

        # FIXME this does not work as expected
        self.ui.run_btn.setEnabled(False)

        if os.path.exists(self.snakefile) is False:
            self.warning("snakefile not found. should not happen")
            return

        self.cmd = ['snakemake', "-s", self.snakefile, "--unlock"]
        self.info("Running " + " ".join(self.cmd))
        self.info("Please wait a second. Unlocking working directory")
        # focus on tab with snakemake output
        self.ui.tabs.setCurrentIndex(0)

        self.ui.tabs_pipeline.setEnabled(False)
        try:
            snakemake_proc = sp.Popen(self.cmd, cwd=self.working_dir)
            snakemake_proc.wait()
        except:
            self.critical("Issue while unlocking the directory")
        finally:
            self.ui.tabs_pipeline.setEnabled(True)

        self.info("unlocking done")
        self.output.append('<font style="color:brown">Unlocking working directory</font>')

        self.ui.run_btn.setEnabled(True)
        self.ui.stop_btn.setEnabled(False)

    # -----------------------------------------------------------------------
    # DAG footer button
    # -----------------------------------------------------------------------

    def show_dag(self):
        try:
            # This command should work on various platform, just in case
            # we add a try/except
            if easydev.cmd_exists('dot') is False:
                msg = "**dot** command not found. Use 'conda install "
                cmd += "graphviz' to install it."
                self.warning(msg)
                msg = WarningMessage((msg))
                msg.exec_()
                return
        except:
            pass
        finally:
            self.info("Creating DAG image.")

        if self.snakefile is None:
            self.warning("No snakefile")
            return

        # We just need the basename because we will run it in the wkdir
        snakefile = os.path.basename(self.snakefile)
        snakemake_line = ["snakemake", "-s", snakefile]
        snakemake_line += ["--rulegraph"]
        if self.mode == "generic" and self.configfile:
            # make sure to copy the config file
            snakemake_line += ["--configfile"]
            snakemake_line += [os.path.basename(self.generic_factory.configfile)]
        snakemake_line += self.get_until_starting_option()

        # Where to save the SVG (temp directory)
        svg_filename = self._tempdir.path() + os.sep + "test.svg"

        self.info(snakemake_line)
        self.process1.setWorkingDirectory(self.working_dir)
        self.process1.setStandardOutputProcess(self.process2)

        self.process1.start("snakemake", snakemake_line[1:])
        self.process2.start("dot", ["-Tsvg", "-o", svg_filename])

        self.process1.waitForFinished(50000)
        self.process2.waitForFinished(50000)

        if os.path.exists(svg_filename):
            self.diag = SVGDialog(svg_filename)
            self.diag.show()
        else:
            msg = 'Could not create the DAG file.'
            error = str(self.process1.readAllStandardError())
            msg = CriticalMessage(msg, error)
            msg.exec_()
            return

    def open_report(self):

        pref = self.preferences_dialog.ui
        filename = pref.preferences_options_general_htmlpage_value.text()
        if filename == "":
            filename = QW.QFileDialog.getOpenFileNames(self,
               "Select your HTML report", self.working_dir,
                "HTML files (*.html)")[0]
            if len(filename) and os.path.exists(filename[0]):
                filename = filename[0]
            else:
                self.warning("No valid HTML selected and none specified in the preferences.")
                return
        else: # we have a filename hardcoded in the preferences
            if self.working_dir is None:
                self.error("Working directory not set yet")
                return

            filename = self.working_dir + os.sep + filename
            if os.path.exists(filename) is False:
                self.error("%s page does not exist. Check the preferences dialog." % filename)
                return
            else:
                self.info("Reading and openning %s" % filename)

        dialog = self.preferences_dialog.ui # an alias
        url = "file://" + filename

        # The browser executable itself
        browser = dialog.preferences_options_general_browser_value.currentText()

        if browser == "pyqt5":
            self.browser = Browser(url)
            self.browser.show()
        else:
            try:
                easydev.execute("%s %s" % (browser, url), showcmd=False,
                    verbose=False)
            except:
                msg = CriticalMessage("Error browser",
                    "Could not open %s with %s. Is %s on your system or path ?" %
                    (url, browser,browser))
                msg.exec_()
                return

    def create_form_dict(self, layout):
        def _cleaner(value):
            # This is to save the YAML file correctly since the widgets tend to
            # convert None and empty strings as '""' or "''"
            if value in ['None', None, '', '""', "''"]:
                return None
            else:
                # this tries to convert to a list #issue #515 
                try:
                    return eval(value)
                except:
                    return value

        widgets = (layout.itemAt(i).widget() for i in range(layout.count()))
        form_dict = {w.get_name(): _cleaner(w.get_value()) if w.is_option()
                         else self.create_form_dict(w.get_layout())
                         for w in widgets}
        return form_dict

    def clear_form(self):
        self.clear_layout(self.form)

    """def _load_and_merge_config(self):

        self.critical("Entering Load and Merge method")

        config_file = self.working_dir + os.sep + os.path.basename(self.configfile)
        if os.path.isfile(config_file) is False:
            self.critical("load_and_merge: no config file present. Nothing to do")
            return False

        self.warning("A config file is present in the target directory")

        # this should always work but who knows
        try:
            self.critical("load_and_merge: reading config file")
            cfg = snaketools.SequanaConfig(config_file)
            cfg.cleanup() # set all empty strings and %()s to None
            config_dict = cfg.config
        except Exception as err:
            self.warning(err)
            self.critical("Could not interpret the sequana config file")
            return False

        if set(self.config._yaml_code.keys()) == set(config_dict.keys()):
            self.critical("load_and_merge:  config in directory same as input config. You may import it")
            msg = QW.QMessageBox(
                QW.QMessageBox.Question, "Question",
                "A config file already exists in the working directory.\n" +
                "%s.\n" % self.working_dir +
                "Do you want to overwrite it ? (if not, the existing file is imported)?",
                QW.QMessageBox.Yes | QW.QMessageBox.No,
                self, Qt.Dialog | Qt.CustomizeWindowHint)
            # Yes == 16384
            if msg.exec_() != 16384:
                self.critical('replacing config with content found in the directory')
                print(config_dict)
                self.config._yaml_code.update(config_dict)
                self.create_base_form()
                self.fill_until_starting() #self.rule_list)
                return True
            else:
                self.critical('We will overwrite the existing file')
                return False
        else:
            self.critical("The config file that already exists is different. Nothing done")
            return False
        return
    """

    def eventFilter(self, source, event):
        """ Inactivate wheel event of combobox
        """
        if event.type() == QtCore.QEvent.Wheel and source is self.ui.choice_button:
            return True
        return False

    # ---------------------------------------------------
    #  settings and close
    # ---------------------------------------------------

    def read_settings(self):
        self.info("Reading settings")
        settings = QtCore.QSettings("sequana_gui", "mainapp")
        if settings.value("tab_position") is not None:
            index = settings.value("tab_position")
            self.ui.tabs_pipeline.setCurrentIndex(int(index))

        if settings.value("tab_generic_position") is not None:
            index = settings.value("tab_generic_position")
            self.ui.tabs_generic.setCurrentIndex(int(index))

        if settings.value("tab_sequana_position") is not None:
            index = settings.value("tab_sequana_position")
            self.ui.tabs_sequana.setCurrentIndex(int(index))

        if settings.value("tab_sequana_input_position") is not None:
            index = settings.value("tab_sequana_input_position")
            self.ui.tabWidget.setCurrentIndex(int(index))

    def write_settings(self):
        settings = QtCore.QSettings("sequana_gui", "mainapp")

        # tab snakemake output/logger/ipython
        index = self.ui.tabs_pipeline.currentIndex()
        settings.setValue("tab_position", index)

        index = self.ui.tabs_generic.currentIndex()
        settings.setValue("tab_generic_position", index)

        index = self.ui.tabs_sequana.currentIndex()
        settings.setValue("tab_sequana_position", index)

        index = self.ui.tabWidget.currentIndex()
        settings.setValue("tab_sequana_input_position", index)

    def _close(self):
        self.write_settings()
        # end any process running that may be running
        self.click_stop()

        self._tempdir.remove()
        try:self.browser.close()
        except:pass

    def closeEvent(self, event):
        #Close button (red X)
        self._close()

    def close(self):
        # Menu or ctrl+q
        self._close()
        super().close()


class Options(argparse.ArgumentParser):
    def __init__(self, prog="sequana_gui"):
        usage = """Sequanix (part of Sequana project) is a GUI for running Snakefiles

        For Sequana project, you can pre-filled sections as follows:

            sequanix -p quality_control -w analysis -i .

        to prefill the quality_control pipeline to used the local directory to
        search for input files (fastq.gz) and run the analysis in the working
        directory "analysis"

        For Generic snakefiles:

            sequanix -s SNAKEFILE -c CONFIGFILE -w analysis

        will run the snakefile (with its config file) into a working directory.

        """
        description = """"""
        super(Options, self).__init__(usage=usage, prog=prog,
            description=description, formatter_class=easydev.SmartFormatter)
        group = self.add_argument_group("GENERAL")
        group.add_argument("-w", "--working-directory", dest="wkdir",
            help="Set working directory", default=None)
        group.add_argument("-n", "--no-splash", dest="nosplash",
            action="store_true",
            help="No splash screen")

        group = self.add_argument_group("SEQUANA")
        group.add_argument("-p", "--pipeline", dest="pipeline",
            default=None,
            help="A valid sequana pipeline name")

        group_mut = group.add_mutually_exclusive_group()
        group_mut.add_argument("-i", "--input-directory", dest="input_directory",
            default=None,
            help="input directory where to find the input data")
        group_mut.add_argument("-f", "--input-files", dest="input_files",
            default=None, nargs="*",
            help="input files")
        group.add_argument("-C", "--replace-configfile", dest="sequana_configfile",
            default=None,
            help="Replace default sequana config file with local configfile")

        group = self.add_argument_group("GENERIC PIPELINES")
        group.add_argument("-s", "--snakefile", dest="snakefile",
            default=None,
            help="A valid Snakefile")
        group.add_argument("-c", "--configfile", dest="configfile",
            default=None,
            help="optional config file to be used by the Snakefile")
        group.add_argument("-y", "--schema", dest="schemafile",
            default=None,
            help="optional schema file to check the config file")


def main(args=None):

    if args is None:
        args = sys.argv[:]
    user_options = Options()
    options = user_options.parse_args(args[1:])

    signal.signal(signal.SIGINT, sigint_handler)
    app = QW.QApplication(sys.argv)

    filename = sequana_data("drawing.png", "../gui")

    if options.nosplash:
        app.processEvents()
        sequana = SequanaGUI(user_options=options)
        sequana.show()
    else:
        # Show the splash screen for a few seconds
        splash_pix = QtGui.QPixmap(filename)
        splash = QW.QSplashScreen(splash_pix, Qt.WindowStaysOnTopHint)
        splash.setMask(splash_pix.mask())
        splash.show()

        for i in range(0, 100):
            t = time.time()
            while time.time() < t + 0.5/100.:
                app.processEvents()

        app.processEvents()
        sequana = SequanaGUI(user_options=options)
        sequana.show()
        splash.finish(sequana)

    # Make sure the main window is the active one
    sequana.raise_()
    sequana.activateWindow()

    sys.exit(app.exec_())


if __name__ == "__main__":
    main()

