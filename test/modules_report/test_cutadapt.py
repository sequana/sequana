from sequana.modules_report.cutadapt import CutadaptModule
from sequana.utils import config
from sequana import bedtools

from . import test_dir


def test_cutadapt_module(tmpdir):
    directory = tmpdir.mkdir('test_module')
    config.output_dir = str(directory)
    config.sample_name = 'JB409847'
    c = CutadaptModule(f"{test_dir}/data/test_cutadapt_paired.txt", "TEST", "test.html")

def test_output():
    # Used the PCRFree adapters
    filename = f"{test_dir}/data/test_cutadapt_paired.txt"
    mod = CutadaptModule(filename, "sample_name")
    assert mod.jinja['command'].startswith("cutadapt")
    assert mod.jinja['mode'] == 'Paired-end'
    assert mod.jinja['paired_reads1_with_adapters'] == '273'
    assert mod.jinja['paired_reads2_with_adapters'] == '243'
    assert mod.jinja['paired_reads_kept'] == '2,189'

    # atropos results should be identical except for a few differences (e.g.
    # command starts with atropos) This is for atropos 1.0.23 (txt only)
    # Note that here, we used Nextera 
    filename = f"{test_dir}/data/test_atropos_paired.txt"
    mod = CutadaptModule(filename, "sample_name")
    assert mod.jinja['paired_reads1_with_adapters'] == '197'
    assert mod.jinja['paired_reads2_with_adapters'] == '262'
    assert mod.jinja['paired_reads_kept'] == '2,420'


    # singled -end with text
    mod = CutadaptModule(f"{test_dir}/data/test_atropos_se.txt", "sample_name")
    assert mod.jinja['command'].startswith("atropos")


def test_atropos_paired(tmpdir):
    # This is for the new version of cutadapt with version 1.1
    directory = tmpdir.mkdir('test_module')
    config.output_dir = str(directory)
    config.sample_name = 'JB409847'
    c = CutadaptModule(f"{test_dir}/data/test_atropos_pe.json", "TEST", "test.html")
    assert c.jinja['mode'] == 'Paired-end'
    c = CutadaptModule(f"{test_dir}/data/test_atropos_se.json", "TEST", "test.html")
    assert c.jinja['mode'] == 'Singled-end'



