import os
from sequana.scripts import main
from sequana import sequana_data
import pytest
from unittest.mock import patch

prog = "sequana"


def test_main():
    main.main([prog])


def test_version():
    main.main([prog, '--version'])

def test_init(mocker):
    try:
        # py3
        with mocker.patch('builtins.input', return_value="y"):
            try:
                main.main([prog, '--pipeline', "qualitydummy"])
                assert False
            except:
                assert True
    except:
        # py2
        with mocker.patch('__builtin__.input', return_value="y"):
            main.main([prog, '--pipeline', "quality_control"])
        with mocker.patch('__builtin__.input', return_value="y"):
            try:
                main.main([prog, '--pipeline', "qualitydummy"])
                assert False
            except:
                assert True

def test_help():
    try:
        main.main([prog, '--help'])
        assert False
    except SystemExit:
        pass
    else:
        raise Exception

def test_info(mocker):
    def func(*args, **kwargs):
        pass
    with patch('easydev.onweb', func ):
        main.main([prog, '--info', "quality_control"])

def test_show_pipelines():
    main.main([prog, '--show-pipelines'])

def test_mutually_exclusive():
    try:
        main.main([prog, '--pipeline', 'quality_control', '--info', "quality_control"])
        assert False
    except SystemExit:
        assert True

def test_input(tmpdir):
    directory = tmpdir.mkdir("analysis0")
    name = directory.__str__()

    file1 = sequana_data('Hm2_GTGAAA_L005_R1_001.fastq.gz', 'data')
    pathdata, basename = os.path.split(file1)
    file2 = sequana_data('Hm2_GTGAAA_L005_R2_001.fastq.gz', 'data')
    pattern = "Hm2_GTGAAA_L005*"

    main.main([prog, "--pipeline", "quality_control", "--input-directory",
              pathdata, "--input-pattern", pattern, "--no-adapters" , "--force",
              "--working-directory", name])

def test_missing_adapters():
    file1 = sequana_data('Hm2_GTGAAA_L005_R1_001.fastq.gz', 'data')
    pathdata, basename = os.path.split(file1)
    try:
        main.main([prog, "--pipeline", "quality_control", "--input-pattern",
            "Hm*gz", "--input-directory", pathdata])
        assert False
    except:
        assert True

    try:
        main.main([prog, "--pipeline", "quality_control", "--input-pattern", file1,
            "--adapters", "dummy" "--working-directory", name, "--force"])
        assert False
    except:
        assert True

def test_without_cluster_config(tmpdir):
    directory = tmpdir.mkdir("analysis1")
    name = directory.__str__()
    file1 = sequana_data('Hm2_GTGAAA_L005_R1_001.fastq.gz', 'data')
    main.main([prog, "--pipeline", "rnaseq", "--input-pattern", file1, "--snakemake-cluster", 
        '"sbatch --mem={cluster.ram}"', "--ignore-cluster-config", "--force",
        "--working-directory", name, "--no-adapters"])
    with open("%s/runme.sh" %  name) as fh: 
        assert "cluster_config" not in fh.read()

# Fails on travis
def _test_with_cluster_config(tmpdir):
    directory = tmpdir.mkdir("analysis2")
    name = directory.__str__()
    file1 = sequana_data('Hm2_GTGAAA_L005_R1_001.fastq.gz', 'data')
    main.main([prog, "--pipeline", "rnaseq", "--input-pattern", file1, "--snakemake-cluster", 
        '"sbatch --mem={cluster.ram}"',"--force",
        "--working-directory", name])
    with open("%s/runme.sh" %  name) as fh: 
        assert "cluster_config" in fh.read()


def test_get_config():

    main.main([prog,  "--pipeline", "quality_control",
        "--get-config"])
    import os
    os.remove("config.yaml")



def test_config_params(tmpdir):
    directory = tmpdir.mkdir("analysis3")
    name = directory.__str__()
    file1 = sequana_data('Hm2_GTGAAA_L005_R1_001.fastq.gz', 'data')
    main.main([prog, "--pipeline", "quality_control", "--force", "--input-pattern", file1,
        "--no-adapters",  "--config-params",  "bwa_mem_phix:threads:4" ,
        "--working-directory", name])

    try:
        main.main([prog, "--pipeline", "quality_control", "--force",  "--pattern", file1,
            "--no-adapters",  "--config-params",  "bwa_mem_phix",
            "--working-directory", name])
        assert False
    except:
        assert True


def test_adapters(tmpdir):
    directory = tmpdir.mkdir("analysis3")
    name = directory.__str__()
    file1 = sequana_data('Hm2_GTGAAA_L005_R1_001.fastq.gz', 'data')
    main.main([prog, "--pipeline", "quality_control", "--input-pattern", file1,
        "--adapters", "Nextera", "--working-directory", name, "--force"])

    main.main([prog, "--pipeline", "quality_control", "--input-pattern", file1,
        "--adapter-fwd", "ACGTACGT", "--adapter-rev", "ACGTACGT",
        "--working-directory", name, "--force"])

    main.main([prog, "--pipeline", "quality_control", "--input-pattern", file1,
        "--adapters", "universal", "--working-directory", name, "--force"])








