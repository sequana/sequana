from sequana import snaketools, sequana_data
from sequana.snaketools import DOTParser
import os, shutil
import tempfile
from sequana import Module, SequanaConfig


def test_dot_parser():
    s = DOTParser(sequana_data("test_dag.dot", "testing"))
    s.add_urls()
    try:os.remove("test_dag.ann.dot")
    except:pass

def test_md5():
    from sequana import Module
    m = Module("quality_control")
    data = m.md5()

def test_modules():
    assert "dag" in snaketools.modules.keys()
    assert snaketools.modules['dag'].endswith("dag.rules")


def test_getcleanup_rules():
    filename =  snaketools.modules['fastq_sampling']
    try:
        snaketools.get_cleanup_rules(filename)
    except:
        pass


def test_snakemake_stats():
    # this is created using snakemake with the option "--stats stats.txt"
    s = snaketools.SnakeMakeStats(sequana_data("test_snakemake_stats.txt"))
    s.plot()


def test_module():
    # a rule without README
    m = snaketools.Module('mark_duplicates')
    m.description
    print(m)
    m.path
    m.snakefile

    # a rule with README
    m = snaketools.Module('dag')
    m.description
    m.overview
    m.is_executable()
    m.check()

    # a pipeline
    m = snaketools.Module('quality_control')
    m.is_executable()
    m.check()
    m.snakefile
    m.name
    m
    print(m)
    assert m.cluster_config.endswith("cluster_config.json")
    assert m.schema_config.endswith("schema.yaml")


def _test_module_onweb():
    m = snaketools.Module('quality_control')
    m.onweb()


def test_valid_config():
    config = snaketools.SequanaConfig(None)

    s = snaketools.Module("quality_control")
    config = snaketools.SequanaConfig(s.config)

    from easydev import TempFile
    with TempFile() as fh:
        config.save(fh.name)


def test_sequana_config():
    s = snaketools.Module("quality_control")
    config = snaketools.SequanaConfig(s.config)

    assert config.config.get("kraken:dummy", "test") == "test"
    assert config.config.get("kraken:dummy") == None

    # --------------------------------- tests different constructors
    config = snaketools.SequanaConfig()
    config = snaketools.SequanaConfig({"test":1})
    assert config.config.test == 1
    # with a dictionary
    config = snaketools.SequanaConfig(config.config)
    # with a sequanaConfig instance
    config = snaketools.SequanaConfig(config)
    # with a non-yaml file
    try:
        json = sequana_data('test_summary_fastq_stats.json')
        config = snaketools.SequanaConfig(json)
        assert False
    except:
        assert True
    try:
        config = snaketools.SequanaConfig("dummy_dummy")
        assert False
    except:
        assert True

    # Test an exception
    s = snaketools.Module("quality_control")
    config = snaketools.SequanaConfig(s.config)
    config._recursive_update(config._yaml_code, {"input_directory_dummy": "test"})

    # loop over all pipelines, read the config, save it and check the content is
    # identical. This requires to remove the templates. We want to make sure the
    # empty strings are kept and that "no value" are kept as well
    #
    #    field1: ""
    #    field2:
    #
    # is unchanged
    from easydev import TempFile
    output = TempFile(suffix=".yaml")
    for pipeline in snaketools.pipeline_names:
        config_filename = Module(pipeline)._get_config()
        cfg1 = SequanaConfig(config_filename)
        cfg1.cleanup() # remove templates and strip strings

        cfg1.save(output.name)
        cfg2 = SequanaConfig(output.name)
        assert cfg2._yaml_code == cfg1._yaml_code
        cfg2._update_config()
        assert cfg1.config == cfg2.config
    output.delete()


def test_message():
    snaketools.message("test")


def test_dummy_manager():
    ss = snaketools.DummyManager()
    ss = snaketools.DummyManager(["test1.fastq.gz", "test2.fastq.gz"])
    assert ss.paired == True
    ss = snaketools.DummyManager(["test1.fastq.gz"])
    assert ss.paired == False
    ss = snaketools.DummyManager("test1.fastq.gz")
    assert ss.paired == False


def test_pipeline_manager():

    # test missing input_directory
    cfg = SequanaConfig({})
    try:
        pm = snaketools.PipelineManager("custom", cfg)
        assert False
    except:
        assert True

    # normal behaviour but no input provided:
    config = Module("quality_control")._get_config()
    cfg = SequanaConfig(config)
    cfg.cleanup() # remove templates
    try:
        pm = snaketools.PipelineManager("custome", cfg)
        assert False
    except:
        assert True

    cfg = SequanaConfig(config)
    cfg.cleanup() # remove templates
    file1 = sequana_data("Hm2_GTGAAA_L005_R1_001.fastq.gz")
    cfg.config.input_directory, cfg.config.input_pattern = os.path.split(file1)
    #file2 = sequana_data("Hm2_GTGAAA_L005_R2_001.fastq.gz")
    pm = snaketools.PipelineManager("custom", cfg)
    assert pm.paired == False

    cfg = SequanaConfig(config)
    cfg.cleanup() # remove templates
    cfg.config.input_directory, cfg.config.input_pattern = os.path.split(file1)
    cfg.config.input_pattern = "Hm*gz"
    #file1 = sequana_data("Hm2_GTGAAA_L005_R1_001.fastq.gz")
    pm = snaketools.PipelineManager("custom", cfg)
    assert pm.paired == True

    pm.getlogdir("fastqc")
    pm.getwkdir("fastqc")
    pm.getrawdata()
    pm.getreportdir("test")
    pm.getname("fastqc")


def test_file_name_factory():
    import glob

    def inner_test(ff):
        len(ff)
        print(ff)
        ff.filenames
        ff.realpaths
        ff.all_extensions
        ff.pathnames
        ff.extensions

    #list
    list_files = glob.glob("*.py")
    ff = snaketools.FileFactory(list_files)
    inner_test(ff)

    # glob
    ff = snaketools.FileFactory("*py")
    inner_test(ff)


    directory = os.path.dirname(sequana_data("Hm2_GTGAAA_L005_R1_001.fastq.gz"))

    ff = snaketools.FastQFactory(directory + "/Hm2*fastq.gz", verbose=True)
    assert ff.tags == ['Hm2_GTGAAA_L005']

    ff.get_file1(ff.tags[0])
    ff.get_file2(ff.tags[0])
    assert len(ff) == 1


def test_copy_requirements():
    # We need 4 cases:
    # 1- http 
    # 2- a sequana file (phix)
    # 3- an existing file elsewhere (here just a temporary file)
    # 4- an existing file in the same directory as the target dir

    from easydev import TempFile
    fh = tempfile.TemporaryDirectory()
    targetdir = fh.name

    # Case 3: a temporary file
    temprequire = TempFile()

    # Case 4: a local file (copy of the temp file) 
    # TODO
    #localfile = temprequire.name.split(os.sep)[-1]
    #shutil.copy(temprequire.name, targetdir)

    cfg = snaketools.SequanaConfig()
    cfg.config.requirements = ["phiX174.fa", temprequire.name, 
        #localfile,
        "https://raw.githubusercontent.com/sequana/sequana/master/README.rst"]
    cfg._update_yaml()
    cfg.copy_requirements(target=fh.name)

    # error
    cfg.config.requirements = ['dummy']
    try:
        cfg.copy_requirements(target=fh.name)
        assert False
    except:
        assert True


def test_onsuccess(tmpdir):
    directory = tmpdir.mkdir("onsuccess")
    p1 = directory.join("Makefile")
    p2 = directory.join("cleanup.py")
    onsuc = snaketools.OnSuccess()
    onsuc.makefile_filename = p1
    onsuc.makefile_cleanup = p2


def test_onsuccess_cleaner():
    fh = tempfile.TemporaryDirectory()
    onsucc = snaketools.OnSuccessCleaner()
    onsucc.makefile_filename = fh.name + os.sep + "Makefile"
    onsucc.add_bundle()
    onsucc.add_makefile()


def test_build_dynamic_rule():

    code = "whatever"
    fh = tempfile.TemporaryDirectory()
    directory = fh.name
    snaketools.build_dynamic_rule(code, directory)


def test_init():
    snaketools.init("quality_control.rules", globals())
    assert "expected_output" in globals()





def test_get_pipeline_statistics():
    df = snaketools.get_pipeline_statistics()

