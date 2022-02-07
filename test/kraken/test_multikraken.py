from sequana import MultiKrakenResults2, MultiKrakenResults

from . import test_dir


def test_mkr(tmpdir):
    mkr = MultiKrakenResults(
        [
            f"{test_dir}/data/test_kraken_multiple_1.csv",
            f"{test_dir}/data/test_kraken_multiple_2.csv",
        ],
        sample_names=["a", "b"],
    )

    mkr = MultiKrakenResults(
        [
            f"{test_dir}/data/test_kraken_multiple_1.csv",
            f"{test_dir}/data/test_kraken_multiple_2.csv",
        ]
    )

    filename = tmpdir.mkdir("temp").join("test.png")
    mkr.plot_stacked_hist(kind="bar", output_filename=filename)
    mkr.plot_stacked_hist(kind="barh", output_filename=filename)


def test_mkr2(tmpdir):
    mkr = MultiKrakenResults2(
        [
            f"{test_dir}/data/test_kraken_mkr2_summary_1.json",
            f"{test_dir}/data/test_kraken_mkr2_summary_2.json",
        ],
        sample_names=["a", "b"],
    )

    mkr = MultiKrakenResults2(
        [
            f"{test_dir}/data/test_kraken_mkr2_summary_1.json",
            f"{test_dir}/data/test_kraken_mkr2_summary_2.json",
        ]
    )

    filename = tmpdir.mkdir("temp").join("test.png")
    mkr.plot_stacked_hist(output_filename=filename)
