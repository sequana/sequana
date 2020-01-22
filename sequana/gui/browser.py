# coding: utf-8
from PyQt5 import QtCore, QtGui, Qt, QtWidgets

from packaging import version

msg = """WARNING:: Qt import failed. This is a know issue. Qt import should work for version
5.9 and above, and qt 5.6 Other version may not work. If you get this message,
please post an issue on https://github.com/sequana/sequana with the python
version, the qt version and this error message:\n\n"""

if version.parse(QtCore.QT_VERSION_STR) > version.parse("5.9"):
    # In version 5.9, QtWebEngineWidgets and QtWebKit are deprecated
    #from PyQt5.Qt import QTabWidget
    #from PyQt5 import QtWebEngine as QtWebKit
    from PyQt5.QtWebEngineWidgets import QWebEngineView as QWebView
    from PyQt5.QtWebEngineWidgets import QWebEnginePage as QWebPage
    from PyQt5.QtWebEngineWidgets import QWebEngineSettings as QWebSettings
    #try:
    #    from PyQt5.Qt import QWebEnginePage as QWebPage
    #    from PyQt5.QtWebKitWidgets import QWebPage
    #except: 
else: # work in PyQt5.6
    try:
        from PyQt5.QtWebKitWidgets import QWebView
        from PyQt5.QtWebKit import QWebSettings
        from PyQt5.Qt import QWebPage
    except Exception as err:
        print(msg + str(err))
        pass



from PyQt5.QtWidgets import QProgressBar, QLineEdit

# potential resources for improvements:
# https://github.com/ralsina/devicenzo/blob/master/devicenzo.py


class Browser(Qt.QMainWindow):
    """

    On purpose, there is no caching so that (if re-generated), the
    new content of an HTML is shown.

    """
    def __init__(self, url):
        Qt.QMainWindow.__init__(self)

        # Progress bar
        # ------------------------------------------------------------
        self.progress = 0
        # Main page QWebView
        # -------------------------------------------------------------
        self.wb = SequanaQWebView(
            parent=self,
            titleChanged=self.setWindowTitle)
        self.wb.urlChanged.connect(lambda u: self.url.setText(u.toString()))

        self.wb.titleChanged.connect(self.adjustTitle)
        self.wb.loadProgress.connect(self.setProgress)


        self.setCentralWidget(self.wb)

        # Main menu tool bar
        # -------------------------------------------------------------
        self.tb = self.addToolBar("Main Toolbar")
        for a in (  QWebPage.Back,
                    QWebPage.Forward,
                    QWebPage.Reload,
                    QWebPage.DownloadLinkToDisk):
            self.tb.addAction(self.wb.pageAction(a))

        self.url = QLineEdit(returnPressed =lambda:self.wb.setUrl(
            QtCore.QUrl.fromUserInput(self.url.text())))
        self.tb.addWidget(self.url)

        # status bar ---------------------------------------------------
        self.sb = self.statusBar()
        try:
            #pyqt5.6
            self.wb.statusBarMessage.connect(self.sb.showMessage)
        except:
            pass
        self.wb.page().linkHovered.connect(lambda l: self.sb.showMessage(l, 3000))

        # Search bar
        # ------------------------------------------------------------
        self.search = QLineEdit(returnPressed = lambda: self.wb.findText(self.search.text()))
        self.search.show()
        self.search.hide() # To make ctrl+F effective, need to show/hide ?

        # The shortcuts
        # ---------------------------------------------------------
        self.showSearch = Qt.QShortcut("Ctrl+F", self,
            activated= lambda: (self.search.show() , self.search.setFocus()))
        self.hideSearch = Qt.QShortcut("Esc", self,
            activated= lambda: (self.search.hide(), self.wb.setFocus()))
        self.quit = Qt.QShortcut("Ctrl+Q", self, activated=self.close)
        self.zoomIn = Qt.QShortcut("Ctrl++", self,
            activated= lambda: self.wb.setZoomFactor(self.wb.zoomFactor()+.2))
        self.zoomOut = Qt.QShortcut("Ctrl+-", self,
            activated= lambda: self.wb.setZoomFactor(self.wb.zoomFactor()-.2))
        self.zoomOne = Qt.QShortcut("Ctrl+=",
            self, activated = lambda: self.wb.setZoomFactor(1))

        # Add alt+left and right keys to navigate backward and forward
        Qt.QShortcut(QtCore.Qt.AltModifier + QtCore.Qt.Key_Left, self,
            activated=lambda: self.wb.back())
        Qt.QShortcut(QtCore.Qt.AltModifier + QtCore.Qt.Key_Right, self,
            activated=lambda: self.wb.forward())

        # Add components on the page
        self.sb.addPermanentWidget(self.search)

        # Finally, load the URL
        self.wb.load(QtCore.QUrl(url))

        try:self.wb.settings().setObjectCacheCapacities(0,0,0)
        except:pass

    def adjustTitle(self):
        if 0 < self.progress < 100:
            self.setWindowTitle("%s (%s%%)" % (self.wb.title(), self.progress))
        else:
            self.setWindowTitle(self.wb.title())

    def setProgress(self, p):
        self.progress = p
        self.adjustTitle()


class SequanaQWebView(QWebView):
    """This is the webview for the application.

    It represents a browser window, either the main one or a popup.
    It's a simple wrapper around QWebView that configures some basic settings.
    """
    def __init__(self, parent=None, **kwargs):
        """Constructor for the class"""
        super(SequanaQWebView, self).__init__(parent)
        self.kwargs = kwargs

        # Javascript and other settings
        # ------------------------------------------------------------
        try:
            self.settings().setAttribute(
                QWebSettings.JavascriptCanOpenWindows, True)
            self.settings().setAttribute(
                QWebSettings.LocalStorageEnabled, True)
            self.settings().setAttribute(
                QWebSettings.PluginsEnabled, True)
        except:
            print("QtWebKit.QWebSettings not available for you PyQt version")

    def createWindow(self, type):
        """Handle requests for a new browser window.

        Method called whenever the browser requests a new window
        (e.g., <a target='_blank'> or window.open()).
        Overridden from QWebView to allow for popup windows, if enabled.
        """
        #this = Browser(self.url())
        #this.show()

        self.popup = SequanaQWebView(**self.kwargs)
        self.popup.setObjectName("web_content")
        self.popup.setWindowTitle("Sequana browser")
        self.popup.page().windowCloseRequested.connect(self.popup.close)
        self.popup.show()
        return self.popup
