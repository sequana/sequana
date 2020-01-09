from sequana.gui.preferences import PreferencesDialog
import pytest



def test_preference(qtbot, mocker):
    widget = PreferencesDialog()
    qtbot.addWidget(widget)

    widget.read_settings() ##### do not write settings in the test or with a mock

    widget.get_settings()
    widget.accept()
    widget.reject()
