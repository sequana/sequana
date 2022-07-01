from sequana.scripts import start_pipeline
import os
import pytest

prog = "sequana_start_pipeline"


def test_analysis(tmpdir):

    directory = tmpdir.mkdir("temp")

    start_pipeline.main([prog, '--name', "test", "--force", "--no-interaction", 
        "--output-directory", str(directory), "--description", "whatever", "--keywords", "1,2,3"])


def test_help():
    with pytest.raises(SystemExit):
        start_pipeline.main([prog, '--help', '1>/tmp/out', '2>/tmp/err'])

