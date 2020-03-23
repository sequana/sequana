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
import colorlog
import shutil

from PyQt5 import QtWidgets as QW
from PyQt5 import QtCore

__all__ = ['Logger', 'Tools', 'QPlainTextEditLogger']


class Logger(object):
    """Aliases to colorlog different methods (e.g. info, debug)

    Set a stream handler to the filename set in the constructor, which 
    can be changed using the attribute :attr:`_logger_output`

    """
    def __init__(self, filename="sequana_logger_debug.txt"):
        self._logger_output = filename
        self.init_logger()

    def init_logger(self):
        self._mylogger = colorlog.getLogger("sequanix")
        """self._fh = open(self._logger_output, "w")
        self._handler = colorlog.StreamHandler(self._fh)
        formatter = colorlog.ColoredFormatter(
            "%(log_color)s%(levelname)-8s%(reset)s %(blue)s%(message)s",
            datefmt=None,
            reset=True,
            log_colors={
                'DEBUG':    'cyan',
                'INFO':     'green',
                'WARNING':  'yellow',
                'ERROR':    'red',
                'CRITICAL': 'red',
            }
        )
        self._handler.setFormatter(formatter)
        self._mylogger.addHandler(self._handler)
        """
    def save_logger(self):
        self._handler.close()
        self._fh.close()
        return self._logger_output
        #self.init_logger()

    def info(self, text):
        self._mylogger.info(text)

    def error(self, text):
        self._mylogger.error(text)

    def debug(self, text):
        self._mylogger.debug(text)

    def critical(self, text):
        self._mylogger.critical(text)

    def warning(self, text):
        self._mylogger.warning(text)


class Tools(Logger):
    def __init__(self):
        super(Tools, self).__init__()

    def copy(self, source, target):
        try:
            shutil.copy(source, target)
        except Exception as err:
            self.error(err)
            self.warning("Cannot overwrite existing file. (Probably identical)")


class QPlainTextEditLogger(colorlog.StreamHandler, QtCore.QObject):
    """


    March 2020. On travis and locally, some tests failed randomly
    with a RuntimeError: wrapped c/c++ object of type qplaintextedit has been deleted

    IT happens to be a issue in this class. According to this stackoverflow
    entry, https://stackoverflow.com/questions/28655198/best-way-to-display-logs-in-pyqt
    we had to implement a thread-safe version for travis and pytest to run
    properly. 

    """
    appendHtml = QtCore.pyqtSignal(str)

    def __init__(self, parent):
        super().__init__()
        self.widget = QW.QPlainTextEdit(parent)
        QtCore.QObject.__init__(self)

        self.widget.setReadOnly(True)
        self.bgcolor = "#aabbcc"
        self.widget.setStyleSheet("background-color: %s" % self.bgcolor)

        self.widget.setStyleSheet("""* {
            selection-background-color: #5964FF;
            background-color: %s
            }""" % self.bgcolor);

        self.appendHtml.connect(self.widget.appendHtml)

    def emit(self, record):
        formatter = """<span style="color:%(color)s;
                        font-weight:%(weight)s">%(msg)s</span>"""
        # "\e[1;31m This is red text \e[0m"
        self.record = record
        msg = self.format(record)
        msg = msg.rstrip("\x1b[0m")
        if msg.startswith('\x1b[31m\x1b[47m'): # critical
            msg = msg.replace("\x1b[31m\x1b[47m", "")
            params = {'msg':msg, 'weight':"bold", "color":"red"}
            #self.widget.appendHtml(formatter % params)
        elif msg.startswith('\x1b[32m'): # info
            msg = msg.replace("\x1b[32m", "")
            params = {'msg':msg, 'weight':"normal", "color":"green"}
            #self.widget.appendHtml(formatter % params)
        elif msg.startswith('\x1b[33m'): # warning
            msg = msg.replace("\x1b[33m", "")
            params = {'msg':msg, 'weight':"normal", "color":"yellow"}
            #self.widget.appendHtml(formatter % params)
        elif msg.startswith('\x1b[31m') or msg.startswith("\\x1b[31m"): # error
            msg = msg.replace("\x1b[31m", "")
            msg = msg.replace("\\x1b[31m", "")
            params = {'msg':msg, 'weight':"normal", "color":"red"}
            #self.widget.appendHtml(formatter % params)
        elif msg.startswith('\x1b[36m'): # debug
            msg = msg.replace("\x1b[36m", "")
            params = {'msg':msg, 'weight':"normal", "color":"cyan"}
            #self.widget.appendHtml(formatter % params)
        else:
            pass
            #self.widget.appendHtml(msg)
        self.appendHtml.emit(formatter % params)
        self.msg = msg

