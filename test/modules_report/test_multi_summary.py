from sequana.modules_report import multi_summary


def test_multi_summary():

    ms = multi_summary.MultiSummary("summary_qc.json")
    ms.get_trimming_percent()
    ms.create_report_content()
