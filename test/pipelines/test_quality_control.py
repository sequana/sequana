import os
import tempfile
from sequana import sequana_data
from easydev import execute


def test_running():
    data = sequana_data("Hm2_GTGAAA_L005_R1_001.fastq.gz")
    input_directory = os.path.dirname(data)

    with tempfile.TemporaryDirectory() as output_directory:
        cmd = "sequana --pipeline quality_control --input-directory %s " + \
            "--output-directory %s --no-adapters --force" 
        execute(cmd % (input_directory, output_directory))
        print(output_directory)
        import time
        time.sleep(10)

        cwd = os.getcwd()
        os.chdir(output_directory)
        execute("sh runme.sh")
        # FIXME: more functional test. This is enough for now to check that the
        # pipeline ends correctly
        assert os.path.exists('Hm2_GTGAAA_L005/report_qc_Hm2_GTGAAA_L005/summary.json')
        os.chdir(cwd)


test_running()
