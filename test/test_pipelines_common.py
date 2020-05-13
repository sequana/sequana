from sequana.pipelines_common import CutadaptOptions
from easydev import AttrDict
import argparse

from sequana import sequana_data




def test_cutadapt_options():

    p = argparse.ArgumentParser()
    so = CutadaptOptions()
    so.add_options(p)

    # test the adapter choice
    for this in ["universal", "PCRFree", "none"]:
        options = {
            "cutadapt_adapter_choice": this,
            "cutadapt_design_file": None,
            "cutadapt_fwd": None,
            "cutadapt_rev": None,
            "skip_cutadapt": False,
        }
        #p.parse_args([])
        options = AttrDict(**options)
        so.check_options(options)


    # test for a valid design and adapter choice
    options = {
        "cutadapt_adapter_choice": "TruSeq",
        "cutadapt_design_file": sequana_data("test_expdesign_Hm2.csv"),
        "cutadapt_fwd": None,
        "cutadapt_rev": None,
        "skip_cutadapt": False,
    }
    options = AttrDict(**options)
    so.check_options(options)

    # test for a valid design but wrong adapter choice
    options = {
        "cutadapt_adapter_choice": "Nextera",
        "cutadapt_design_file": sequana_data("test_expdesign_Hm2.csv"),
        "cutadapt_fwd": None,
        "cutadapt_rev": None,
        "skip_cutadapt": False,
    }
    options = AttrDict(**options)
    try:
        so.check_options(options)
        assert False
    except:
        assert True

    # wrong combo (missing adapter choice)
    options = {
        "cutadapt_adapter_choice": None,
        "cutadapt_design_file": sequana_data("test_expdesign_Hm2.csv"),
        "cutadapt_fwd": None,
        "cutadapt_rev": None,
        "skip_cutadapt": False,
    }
    options = AttrDict(**options)
    try:
        so.check_options(options)
        assert False
    except: assert True

    # wrong quality (missing adapter choice)
    try:
        p.parse_args(["--cutadapt-quality", "-1"])
        assert False
    except: assert True
    p.parse_args(["--cutadapt-quality", "10"])

    # test for a valid design and adapter choice but also fwd/rev provided
    # whereas, we cannot do anything with this combo
    options = {
        "cutadapt_adapter_choice": "TruSeq",
        "cutadapt_design_file": sequana_data("test_expdesign_Hm2.csv"),
        "cutadapt_fwd": "ACGT",  # dummy values
        "cutadapt_rev": "CGTA",  # dummy values
        "skip_cutadapt": False,
    }
    options = AttrDict(**options)
    try:
        so.check_options(options)
        assert False
    except: assert True

    options = {
        "cutadapt_adapter_choice": None,
        "cutadapt_design_file": None,
        "cutadapt_fwd": sequana_data("TruSeqCD_DNA_fwd.fa"), 
        "cutadapt_rev": sequana_data("TruSeqCD_DNA_rev.fa"),
        "skip_cutadapt": False,
    }
    options = AttrDict(**options)
    so.check_options(options)


def test_snakemake_options():
    from sequana.pipelines_common import SnakemakeOptions
    p = argparse.ArgumentParser()
    so = SnakemakeOptions()
    so.add_options(p)
    p.parse_args([])


def test_krakenl_options():
    from sequana.pipelines_common import KrakenOptions
    p = argparse.ArgumentParser()
    so = KrakenOptions()
    so.add_options(p)
    p.parse_args([])

def test_slurm_options():
    from sequana.pipelines_common import SlurmOptions
    p = argparse.ArgumentParser()
    so = SlurmOptions()
    so.add_options(p)
    p.parse_args([])

def test_input_options():
    from sequana.pipelines_common import InputOptions
    p = argparse.ArgumentParser()
    so = InputOptions()
    so.add_options(p)
    p.parse_args([])


def test_general_options():
    from sequana.pipelines_common import GeneralOptions
    p = argparse.ArgumentParser()
    so = GeneralOptions()
    so.add_options(p)
    p.parse_args([])

def test_print_version():
    from sequana.pipelines_common import print_version
    try:print_version("quality_control")
    except:pass
    try:print_version("sequana_dummy")
    except:pass

def test_get_pipeline_location():
    from sequana.pipelines_common import get_pipeline_location
    get_pipeline_location("quality_control")


def test_pipeline_manager():

    from sequana.pipelines_common import PipelineManager
    from sequana.pipelines_common import SequanaManager
    from sequana.pipelines_common import InputOptions
    from sequana.pipelines_common import SlurmOptions
    from sequana.pipelines_common import SnakemakeOptions
    from sequana.pipelines_common import GeneralOptions
    from sequana.pipelines_common import CutadaptOptions

    p = argparse.ArgumentParser()
    so = GeneralOptions()
    so.add_options(p)

    so = InputOptions()
    so.add_options(p)

    so = SnakemakeOptions()
    so.add_options(p)

    so = SlurmOptions()
    so.add_options(p)

    so = CutadaptOptions()
    so.add_options(p)

    from tempfile import TemporaryDirectory
    with TemporaryDirectory() as fout:
        options = p.parse_args(["--working-directory", fout, "--force"])

        pm = SequanaManager(options, "quality_control")


        from sequana.pipelines_common import get_pipeline_location as getpath
        sharedir = getpath('quality_control')


        pm.config.config.input_directory = sharedir
        pm.config.config.input_readtag = options.input_readtag
        pm.config.config.input_pattern = options.input_pattern
        pm.setup()
        pm.update_config(pm.config.config, options, "cutadapt")
        pm.teardown()


    options = p.parse_args(["--working-directory", fout, "--force", 
        "--run-mode", "local", "--slurm-queue", "common"])
    pm = PipelineManager(options, "quality_control")
    pm.setup()

    options = p.parse_args(["--working-directory", fout, "--force", 
        "--run-mode", "slurm", "--slurm-queue", "common"])
    pm = PipelineManager(options, "quality_control")
    pm.setup()


