from sequana.scripts import lane_merging
from sequana import sequana_data
import os
import pytest

prog = "sequana_lane_merging"


def test_analysis():
    file1 = sequana_data("test_L001_R1_001.fastq.gz")
    dirname = os.path.dirname(file1)

    from tempfile import TemporaryDirectory
    directory = TemporaryDirectory()
    # Test that database must be provided
    try:
        # we try twice to test when --force  is not used
        df = lane_merging.main([prog, '--pattern', dirname + "/test_L00?_R?_001*fastq.gz", 
            "--lanes", "1", "2", "--force", "--output-directory", directory.name ])
        df = lane_merging.main([prog, '--pattern', dirname + "/test_L00?_R?_001*fastq.gz",
            "--lanes", "1", "2", "--force", "--output-directory", directory.name ])
        assert False
    except:
        assert True



    
    df = lane_merging.main([prog, '--pattern', dirname + "/test_L00?_R?_001*fastq.gz", 
        "--lanes", "1", "2", "--force", "--output-directory", directory.name ])



def test_help():
    try:
        lane_merging.main([prog, '--help', '1>/tmp/out', '2>/tmp/err'])
        assert False
    except SystemExit:
        pass
    else:
        raise Exception

