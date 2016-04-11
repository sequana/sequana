from sequana.report_main import SequanaReport


def test_report_base():
    from sequana.report_main import BaseReport    
    br = BaseReport("main/index.html")
    try:
        br.parse()
        assert False
    except:
        assert True

def test_report():

    from sequana import sequana_data
    stats = sequana_data("test_snakemake_stats.txt", "testing")

    from sequana import Module
    module = Module('fix_removal')
    

    r = SequanaReport(snakefile=module.snakefile, configfile=module.config,
        stats=stats)
    r.create_report()
