import sys
import os

from sequana.gui.file_browser import FileBrowser, DirectoryDialog

from sequana import sequana_data
from PyQt5 import QtCore

import pytest
from PyQt5 import QtWidgets as QW


# How to use that ? 
# http://pytest-qt.readthedocs.io/en/latest/tutorial.html
def test_basic_search(qtbot, tmpdir):
    '''
    test to ensure basic find files functionality is working.
    '''
    tmpdir.join('test1_R1.fastq.gz').ensure()
    tmpdir.join('test1_R2.fastq.gz').ensure()

    tmpdir.join('test2_R1.fastq.gz').ensure()
    tmpdir.join('test2_R2.fastq.gz').ensure()


def test_directory_dialog(qtbot, mocker):

    widget = FileBrowser(paired=False, directory=False, file_filter=None)

    qtbot.addWidget(widget)
    widget.show()
    assert widget.isVisible()
    widget.setup_color()

    widget.set_empty_path()
    assert widget.get_filenames() == ""
    assert widget.path_is_setup() == False

    # Now, we open the dialog, which pops up. We need to close it...
    widget = FileBrowser(paired=True, directory=False, file_filter=None)
    qtbot.addWidget(widget)
    #qtbot.mouseClick(widget.btn, QtCore.Qt.LeftButton)

    widget.Nmax = 10
    widget._set_paired_filenames([ "test1.fastq.gz"])
    qtbot.addWidget(widget)

    widget.Nmax = 30
    widget._set_paired_filenames(["test1.fastq.gz", "test2.fastq.gz"])
    qtbot.addWidget(widget)


def test_directory_dialog_2(qtbot, tmpdir):
    widget = DirectoryDialog(None, "test", str(tmpdir), "*gz")
    qtbot.addWidget(widget)


def test_directory_dialog_3(qtbot, tmpdir):
    widget = FileBrowser(paired=True)
    qtbot.addWidget(widget)
    # requires an interaction
    #widget.browse_paired_file()


def test_setfilenames(qtbot, tmpdir):
    #tmpdir.join('AA_R1_.fastq.gz').ensure()
    widget = FileBrowser(paired=True)
    qtbot.addWidget(widget)
    widget.set_filenames("AA_R1_.fastq.gz")
    widget.set_filenames("AA_R2_AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA.fastq.gz")
    


def test_others(qtbot):
    widget = FileBrowser(paired=True)
    qtbot.addWidget(widget)
    widget.set_enable(True)
    widget.set_enable(False)
    widget.clicked_connect(widget.close)
    #text changed error
    #widget.changed_connect(widget.close)
    
