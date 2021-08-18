from sequana.scripts import coverage
from sequana import sequana_data
import pytest

prog = "sequana_coverage"

@pytest.fixture
def coveragefix():
    import os
    # local nosetests execution
    try:os.remove('README')
    except:pass
    try:os.remove('quality.rules')
    except:pass
    try:os.remove('config.yaml')
    except:pass


def test_version():
    try:
        coverage.main([prog, '--version'])
        assert False
    except SystemExit:
        pass
    else:
        raise Exception


def test_help():
    try:
        coverage.main([prog, '--help'])
        assert False
    except SystemExit:
        pass
    else:
        raise Exception

@pytest.mark.xfail(reason="too slow or service may be down")
def test_input(tmpdir):
    import os
    # Download reference in temporary directory so that it is erased if the test
    # fails.
    directory_data = tmpdir.mkdir("datatemp")
    cwd = os.getcwd()
    try:
        os.chdir(directory_data.__str__())
        coverage.main([prog, '--download-reference', "JB409847"])
        os.system("""sed -i s"/>ENA/>JB409847 /" %s/JB409847.fa """ % directory_data.__str__())

        coverage.main([prog, '--download-genbank', "JB409847"])
    except Exception:
        raise Exception
    finally:
        os.chdir(cwd)

    directory_run = tmpdir.mkdir("report")

    filename = sequana_data('JB409847.bed')
    try:
        coverage.main([prog, '-i', filename, "-o", "--output-directory",
                directory_run.__str__(), 
                "--window-median", "3001", "-r", 
                "%s/JB409847.fa" % directory_data.__str__()])
        assert False
    except Exception as err:
        print(err)
        assert True
    print(os.listdir(directory_run.__str__()))
    assert os.path.exists(directory_run.__str__() + os.sep + "multiqc_report.html")
    try:
        coverage.main([prog, '-i', filename, "-o", "--output-directory",
                       directory_run.__str__(),
                       "--window-median", "3001", "-r",
                       "%s/JB409847.fa" % directory_data.__str__(), '-c', '1'])
        assert False
    except Exception as err:
        print(err)
        assert True
    assert os.path.exists(str(directory_run) + os.sep + 'multiqc_report.html')


