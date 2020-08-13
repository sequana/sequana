from sequana import snaketools, sequana_data
from sequana.snaketools import DOTParser
import os, shutil
import tempfile
from sequana import Module, SequanaConfig
from easydev import TempFile
import subprocess


def test_dot_parser():
    s = DOTParser(sequana_data("test_dag.dot", "testing"))
    s.add_urls(mapper={'bwa_fix': "test.html"})
    try:os.remove("test_dag.ann.dot")
    except:pass

    s.mode = "v1"
    s.add_urls(mapper={'bwa_fix': "test.html"})
    try:os.remove("test_dag.ann.dot")
    except:pass


def test_md5():
    from sequana import Module
    m = Module("compressor")
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
    with TempFile() as fout:
        s.plot_and_save(filename=fout.name, outputdir=None)

    with tempfile.TemporaryDirectory() as tempdir:
        s.plot_and_save(filename="test.png", outputdir=tempdir)


def test_plot_stats():

    with tempfile.TemporaryDirectory() as indir:
        shutil.copy(sequana_data("test_snakemake_stats.txt"), 
            indir + "/stats.txt")
        with tempfile.TemporaryDirectory() as outdir:
            snaketools.plot_stats(indir, outdir)
    snaketools.plot_stats("dummy", "dummy")


def test_module():
    # a rule without README
    m = snaketools.Module('mark_duplicates_dynamic')
    m.description
    print(m)
    m # test __repr__
    m.__repr__()
    m.path
    m.snakefile
    m.overview
    # a rule with README
    m = snaketools.Module('dag')
    m.description
    m.overview
    assert m.is_executable()
    m.check()

    # a pipeline
    m = snaketools.Module('compressor')
    m.is_executable()
    m.check()
    m.snakefile
    m.name
    m
    print(m)
    assert m.cluster_config.endswith("cluster_config.json")
    assert m.schema_config.endswith("schema.yaml")


def test_valid_config():
    config = snaketools.SequanaConfig(None)

    s = snaketools.Module("compressor")
    config = snaketools.SequanaConfig(s.config)

    from easydev import TempFile
    with TempFile() as fh:
        config.save(fh.name)


def test_sequana_config():
    s = snaketools.Module("compressor")
    config = snaketools.SequanaConfig(s.config)

    assert config.config.get("compressor")["source"] == "fastq.gz"
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
    s = snaketools.Module("compressor")
    config = snaketools.SequanaConfig(s.config)
    config._recursive_update(config._yaml_code, {"input_directory_dummy": "test"})

    #config.check_config_with_schema(s.schema_config)
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


def test_check_config_with_schema():
    schema = Module("compressor").schema_config
    SequanaConfig(Module("compressor").config).check_config_with_schema(schema) 


def test_module_version():
    Module("snpeff/1.0").version == "1.0"


def test_message():
    snaketools.message("test")



def test_pipeline_manager():

    # test missing input_directory
    cfg = SequanaConfig({})
    try:
        pm = snaketools.PipelineManager("custom", cfg)
        assert False
    except:
        assert True

    # normal behaviour but no input provided:
    config = Module("compressor")._get_config()
    cfg = SequanaConfig(config)
    cfg.cleanup() # remove templates
    try:
        pm = snaketools.PipelineManager("custom", cfg)
        assert False
    except:
        assert True

    # normal behaviour
    cfg = SequanaConfig(config)
    cfg.cleanup() # remove templates
    file1 = sequana_data("Hm2_GTGAAA_L005_R1_001.fastq.gz")
    cfg.config.input_directory, cfg.config.input_pattern = os.path.split(file1)
    pm = snaketools.PipelineManager("custom", cfg)
    assert pm.paired == False

    cfg = SequanaConfig(config)
    cfg.cleanup() # remove templates
    cfg.config.input_directory, cfg.config.input_pattern = os.path.split(file1)
    cfg.config.input_pattern = "Hm*gz"
    #file1 = sequana_data("Hm2_GTGAAA_L005_R1_001.fastq.gz")
    pm = snaketools.PipelineManager("custom", cfg)
    pm.plot_stats()
    assert pm.paired == True

    pm.getlogdir("fastqc")
    pm.getwkdir("fastqc")
    pm.getrawdata()
    pm.getreportdir("test")
    pm.getname("fastqc")

    # Test different configuration of input_directory, input_readtag,
    # input_pattern
    # Test the _R[12]_ paired
    with tempfile.TemporaryDirectory() as tmpdir:
        cfg = SequanaConfig()
        cfgname = tmpdir + "/config.yaml"
        cfg.config.input_pattern = "*fastq.gz"
        cfg.config.input_directory = tmpdir
        cfg.config.input_readtag = "_R[12]_"
        cfg._update_yaml()
        cfg.save(cfgname)
        cmd = "touch {}/test_R1_.fastq.gz".format(tmpdir)
        subprocess.call(cmd.split())
        cmd = "touch {}/test_R2_.fastq.gz".format(tmpdir)
        subprocess.call(cmd.split())
        pm = snaketools.PipelineManager("test", cfgname)
        assert pm.paired == True


    # Test the _[12]_ paired 
    with tempfile.TemporaryDirectory() as tmpdir:
        cfg = SequanaConfig()
        cfgname = tmpdir + "/config.yaml"
        cfg.config.input_pattern = "*fastq.gz"
        cfg.config.input_directory = tmpdir
        cfg.config.input_readtag = "_[12]."
        cfg._update_yaml()
        cfg.save(cfgname)
        cmd = "touch {}/test_1.fastq.gz".format(tmpdir)
        subprocess.call(cmd.split())
        cmd = "touch {}/test_2.fastq.gz".format(tmpdir)
        subprocess.call(cmd.split())
        pm = snaketools.PipelineManager("test", cfgname)
        assert pm.paired is True

    # Test the _R[12]_ single end
    with tempfile.TemporaryDirectory() as tmpdir:
        cfg = SequanaConfig()
        cfgname = tmpdir + "/config.yaml"
        cfg.config.input_pattern = "*fastq.gz"
        cfg.config.input_directory = tmpdir
        cfg.config.input_readtag = "_R[12]_"
        cfg._update_yaml()
        cfg.save(cfgname)
        cmd = "touch {}/test_R1_.fastq.gz".format(tmpdir)
        subprocess.call(cmd.split())
        pm = snaketools.PipelineManager("test", cfgname)
        assert pm.paired is False

    # Test the _R[12]_ single end
    with tempfile.TemporaryDirectory() as tmpdir:
        cfg = SequanaConfig()
        cfgname = tmpdir + "/config.yaml"
        cfg.config.input_pattern = "*fq.gz" # wrong on purpose
        cfg.config.input_directory = tmpdir
        cfg.config.input_readtag = "_R[12]_"
        cfg._update_yaml()
        cfg.save(cfgname)
        cmd = "touch {}/test_R1_.fastq.gz".format(tmpdir)
        subprocess.call(cmd.split())
        try:
            pm = snaketools.PipelineManager("test", cfgname)
            assert False
        except:
            assert True

    # Test the _R[12]_ single end
    with tempfile.TemporaryDirectory() as tmpdir:
        cfg = SequanaConfig()
        cfgname = tmpdir + "/config.yaml"
        cfg.config.input_pattern = "*fastq.gz" 
        cfg.config.input_directory = tmpdir
        cfg.config.input_readtag = "R[12]_"
        cfg._update_yaml()
        cfg.save(cfgname)
        cmd = "touch {}/testR1_.fastq.gz".format(tmpdir)
        subprocess.call(cmd.split())
        cmd = "touch {}/testR2_.fastq.gz".format(tmpdir)
        subprocess.call(cmd.split())
        try:
            pm = snaketools.PipelineManager("test", cfgname)
            assert False
        except:
            assert True


def test_pipeline_manager_generic():
    cfg = SequanaConfig({})
    file1 = sequana_data("Hm2_GTGAAA_L005_R1_001.fastq.gz")
    cfg.config.input_directory, cfg.config.input_pattern = os.path.split(file1)
    cfg.config.input_pattern = "Hm*gz"
    pm = snaketools.PipelineManagerGeneric("quality_control", cfg)
    pm.getlogdir("fastqc")
    pm.getwkdir("fastqc")
    pm.getrawdata()
    pm.getreportdir("test")
    pm.getname("fastqc")
    gg = globals()
    gg['__snakefile__'] = "dummy"
    pm.setup(gg)
    del gg['__snakefile__']
    class WF():
        included_stack = ["dummy", 'dummy']
    wf = WF()
    gg['workflow'] = wf
    pm.setup(gg)
    pm.teardown()

    with tempfile.TemporaryDirectory() as dd:
        multiqc = open(dd + "/multiqc.html", "w")
        multiqc.write("test")
        multiqc.close()
        newfile = dd + "/multiqc.html_tmp_"
        pm.clean_multiqc(dd + "/multiqc.html")

def test_file_name_factory():
    import glob

    def inner_test(ff):
        len(ff)
        print(ff)
        ff.filenames
        ff.realpaths
        ff.all_extensions
        ff.pathnames
        ff.pathname
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
    snaketools.init("compressor.rules", globals())
    assert "expected_output" in globals()


def test_get_pipeline_statistics():
    df = snaketools.get_pipeline_statistics()


def test_create_cleanup():
    with tempfile.TemporaryDirectory() as fout:
        snaketools.create_cleanup(fout)


def test_fastqfactory():
    try:
        snaketools.FastQFactory("*", read_tag='error') 
        assert False
    except:
        assert True

    try:
        snaketools.FastQFactory("*", read_tag='[12]') 
        assert False
    except:
        assert True

    directory = os.path.dirname(sequana_data("Hm2_GTGAAA_L005_R1_001.fastq.gz"))

    ff = snaketools.FastQFactory(directory + os.sep + "Hm2*gz", read_tag='R[12]') 
    assert ff.paired is True
    assert ff.tags == ['Hm2_GTGAAA_L005_']

    ff = snaketools.FastQFactory(directory + os.sep + "Hm2*gz", read_tag=None) 
    assert ff.paired is False
    assert sorted(ff.tags) == sorted(['Hm2_GTGAAA_L005_R2_001', 'Hm2_GTGAAA_L005_R1_001'])


def test_makefile():
    with tempfile.TemporaryDirectory() as fout:
        mk = snaketools.Makefile()
        mk.makefile_filename = fout + "/Makefile"
        mk.add_remove_done()
        mk.add_bundle()
        mk.save()


def test_bundle():
    with tempfile.TemporaryDirectory() as fout:
        os = snaketools.OnSuccess()
        os.makefile_filename = fout + "/Makefile"
        os.cleanup_filename = fout + "/sequana_cleanup.py"
        os.add_makefile()
        os.create_recursive_cleanup()











