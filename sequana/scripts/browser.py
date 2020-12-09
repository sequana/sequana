import sys
from sequana.gui.browser import Browser
from PyQt5 import QtWidgets as QW
from easydev import isurl_reachable
import os


def main(args=None): #pragma: no cover
    app = QW.QApplication(sys.argv)

    if len(sys.argv) == 2:
        argument = sys.argv[1]
        if os.path.exists(argument):
            browser = Browser("file:///" + os.path.abspath(argument))
        elif isurl_reachable(argument):
            if argument.startswith("www"):
                print("http://" + argument)
                browser = Browser("http://" + argument)
            else:
                browser = Browser(argument)
        else:
            print("{} not a local file or reachable URL".format(argument))
            browser = Browser("http://sequana.readthedocs.org/en/master")
    else:
        browser = Browser("http://sequana.readthedocs.org/en/master")

    browser.show()
    sys.exit(app.exec_())

if __name__ == "__main__":
     main()

