from sequana.scripts import coverage
import pytest
import os

prog = "sequana_coverage"

from .. import test_dir


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

#@pytest.mark.xfail(reason="too slow or service may be down")
def test_download_reference(tmpdir):
    # Download reference in temporary directory so that it is erased if the test
    # fails.
    directory_data = tmpdir.mkdir("datatemp")
    cwd = os.getcwd()
    os.chdir(directory_data.__str__())
    coverage.main([prog, '--download-reference', "JB409847"])
    os.system("""sed -i s"/>ENA/>JB409847 /" %s/JB409847.fa """ % directory_data.__str__())

    coverage.main([prog, '--download-genbank', "JB409847"])

def test_run(tmpdir):

    directory_run = tmpdir.mkdir("report")

    bedfile = f"{test_dir}/data/bed/JB409847.bed"
    genbank = f"{test_dir}/data/genbank/JB409847.gbk"
    fastafile = f"{test_dir}/data/fasta/JB409847.fasta"

    coverage.main([prog, '-i', bedfile, "-o", "--output-directory",
                directory_run.__str__(), "-b", genbank,
                "--window-median", "3001", "-r",
                fastafile])
    print(os.listdir(directory_run.__str__()))
    #assert os.path.exists(str(directory_run / "multiqc_report.html"))




    coverage.main([prog, '-i', bedfile, "-o", "--output-directory",
                   directory_run.__str__(), "--no-multiqc","--no-html",
                   "--window-median", "3001", "-r",
                   fastafile ])


