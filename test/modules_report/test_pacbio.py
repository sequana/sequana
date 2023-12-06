from easydev import TempFile

from sequana.modules_report import pacbio_input_bam

from . import test_dir


def test_pacbio_input_bam(tmpdir):
    # we need a summary and a bunch of images
    filename = f"{test_dir}/data/summary_pacbio_qc1.json"

    # mock the PNG files found in the summary
    import json

    summary = json.load(open(filename))
    pngname = f"{test_dir}/data/no_data.jpg"
    summary["images"]["gc_vs_length"] = pngname
    summary["images"]["hist_gc_content"] = pngname
    summary["images"]["hist_read_length"] = pngname
    summary["images"]["hist_snr"] = pngname
    summary["images"]["hist_zmw"] = pngname

    summary_file = TempFile()
    with open(summary_file.name, "w") as ff:
        json.dump(summary, ff)

    # Now that we have this new summary file, let us use it
    # we also need an output handler
    ff = TempFile()

    from sequana.utils import config

    config.output_dir = "/tmp"
    # here, ff.name is of the form /tmp/djhfjh4dz so we need to remove the /tmp
    pacbio_input_bam.PacbioInputBAMModule(summary_file.name, ff.name.split("/")[1])

    # cleanup
    summary_file.delete()
    ff.delete()
