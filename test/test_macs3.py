from sequana.macs3 import MACS3Reader, PeakConsensus

from . import test_dir


def test_macs3reader():

    # reading a narrow peak and the .xls file that must be identical
    df1 = MACS3Reader(f"{test_dir}/data/macs3/example_rep1_peaks.narrowPeak")
    df1_bis = MACS3Reader(f"{test_dir}/data/macs3/example_rep1_narrow_peaks.xls")
    assert len(df1) == 5

    # reading a broad peak file
    df2 = MACS3Reader(f"{test_dir}/data/macs3/example_rep1_peaks.broadPeak")

    df1.plot_volcano()
    df1.plot_volcano(plotly=True)


def test_peakconsensus(tmpdir):

    pc = PeakConsensus(f"{test_dir}/data/macs3/example_rep1_peaks.narrowPeak",
                        f"{test_dir}/data/macs3/example_rep2_peaks.narrowPeak")
    df = pc.merge()
    assert len(df) ==7
    pc.plot_venn()


    fout = tmpdir.join("out.saf")
    pc.to_saf(fout)
    fout = tmpdir.join("out.bed")
    pc.to_bed(fout)
