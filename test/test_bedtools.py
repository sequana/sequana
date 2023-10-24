import pytest
from easydev import TempFile

from sequana import bedtools

from . import test_dir


def test_threshold():
    t = bedtools.DoubleThresholds(-5, 5)
    assert t.low == -5
    assert t.high == 5
    assert t.low2 == -2.5
    t = bedtools.DoubleThresholds(-4, 3)
    assert t.low == -4
    assert t.high == 3

    t = bedtools.DoubleThresholds(-8, 8)
    t.ldtr = 0.25
    t.hdtr = 0.25
    assert t.low2 == -2
    assert t.high2 == 2
    print(t)

    t.ldtr = 0.5
    t.hdtr = 0.5
    t.low = -3
    t.high = 3
    assert t.low2 == -1.5
    assert t.high2 == 1.5

    with pytest.raises(Exception):
        t = bedtools.DoubleThresholds(3, 4)

    with pytest.raises(Exception):
        t = bedtools.DoubleThresholds(3, -4)


def test_genomecov():
    filename = f"{test_dir}/data/bed/JB409847.bed"

    # wrong file
    with pytest.raises(Exception):
        bed = bedtools.GenomeCov("dummy.csv")

    # wrong threshold
    with pytest.raises(Exception):
        bed = bedtools.GenomeCov(filename, high_threshold=2)

    # wrong threshold
    with pytest.raises(Exception):
        bed = bedtools.GenomeCov(filename, low_threshold=-2)

    # wrong genbank
    with pytest.raises(SystemExit):
        bed = bedtools.GenomeCov(filename, "dummy.gbk")

    # !now let us read the good data sets by chunkd
    bed = bedtools.GenomeCov(filename, f"{test_dir}/data/genbank/JB409847.gbk", chunksize=5000)
    for c in bed.chr_list:
        c.run(1001, k=2)

    # setter must be bool
    with pytest.raises(SystemExit):
        bed.circular = 1

    # cant use setter
    with pytest.raises(SystemExit):
        bed.feature_dict = {}

    assert len(bed) == 1
    # a getter for the first chromosome
    bed[0]

    # setter available but not sure this is useful
    bed.window_size = 4000
    bed.window_size = 4001
    bed.hist()

    # This requires to call other method before
    for chrom in bed:
        chrom.moving_average(n=501)
        chrom.running_median(n=501, circular=True)
        chrom.running_median(n=501, circular=False)

        chrom.compute_zscore()
        roi = chrom.get_rois()
        with TempFile(suffix=".png") as fh:
            chrom.plot_coverage(filename=fh.name)
        with TempFile(suffix=".png") as fh:
            chrom.plot_hist_zscore(filename=fh.name)
        with TempFile(suffix=".png") as fh:
            chrom.plot_hist_normalized_coverage(filename=fh.name)

        len(chrom)
        print(chrom)
        chrom.get_size()
        chrom.DOC
        chrom.CV
    with TempFile(suffix=".csv") as fh:
        bed.gc_window_size = 100
        bed.to_csv(fh.name)

    # try wrong length:

    with pytest.raises(Exception):
        for chrom in bed:
            chrom.moving_average(n=5000000000000000)

    # plotting
    bed.chr_list[0].plot_hist_coverage()
    bed.chr_list[0].plot_hist_coverage(logx=False, logy=True)
    bed.chr_list[0].plot_hist_coverage(logx=True, logy=False)
    with TempFile(suffix=".png") as fh:
        bed.chr_list[0].plot_hist_coverage(logx=False, logy=False, filename=fh.name)


def test_chromosome():
    filename = f"{test_dir}/data/bed/JB409847.bed"
    # using chunksize of 7000, we test odd number
    bed = bedtools.GenomeCov(filename, f"{test_dir}/data/genbank/JB409847.gbk", chunksize=7000)
    chrom = bed.chr_list[0]
    chrom.run(501, k=2, circular=True)
    print(chrom)

    # using chunksize of 7000, we test even number
    bed = bedtools.GenomeCov(filename, f"{test_dir}/data/genbank/JB409847.gbk", chunksize=7000)
    chrom = bed.chr_list[0]
    chrom.run(501, k=2, circular=True)

    # no chunksize
    bed = bedtools.GenomeCov(filename, f"{test_dir}/data/genbank/JB409847.gbk")
    chrom = bed.chr_list[0]
    chrom.run(501, k=2, circular=True)
    print(chrom)

    # no chunksize
    bed = bedtools.GenomeCov(filename, f"{test_dir}/data/genbank/JB409847.gbk")
    chrom = bed.chr_list[0]
    try:
        chrom._coverage_scaling()
        assert False
    except KeyError:
        assert True
    except:
        assert False

    # zscore not computed yet, so error
    try:
        chrom.plot_rois(3000, 8000)
        assert False
    except:
        assert True
    chrom.run(4001)
    chrom.plot_rois(3000, 8000)


def test_gc_content():
    bed = f"{test_dir}/data/bed/JB409847.bed"
    fasta = f"{test_dir}/data/fasta/JB409847.fasta"
    cov = bedtools.GenomeCov(bed)
    #cov.chrom_names.append("dummy")
    cov.compute_gc_content(fasta)
    cov.get_stats()
    ch = cov[0]
    ch.moving_average(4001, circular=True)
    ch.running_median(4001, circular=True)
    ch.compute_zscore()

    ch.evenness
    ch.CV
    assert ch.get_centralness() > 0.84 and ch.get_centralness() < 0.85
    with TempFile(suffix=".png") as fh:
        ch.plot_gc_vs_coverage(filename=fh.name)

    with TempFile() as fh:
        ch.to_csv(fh.name)

    ch.get_max_gc_correlation(fasta)


def test_ChromosomeCovMultiChunk():
    filename = f"{test_dir}/data/bed/JB409847.bed"
    # using chunksize of 7000, we test odd number
    bed = bedtools.GenomeCov(filename, f"{test_dir}/data/genbank/JB409847.gbk", chunksize=7000)
    bed.compute_gc_content(f"{test_dir}/data/fasta/JB409847.fasta")
    chrom = bed.chr_list[0]
    res = chrom.run(501, k=2, circular=True)
    res.get_summary()
    res.get_rois()


def test_binning():
    filename = f"{test_dir}/data/bed/JB409847.bed"
    # using chunksize of 7000, we test odd number
    bed = bedtools.GenomeCov(filename, f"{test_dir}/data/genbank/JB409847.gbk", chunksize=7000)
    chrom = bed.chr_list[0]
    chrom.run(501, k=2, circular=True, binning=2, cnv_delta=100)


def test_chromosome_intname():
    filename = f"{test_dir}/data/bed/unicycler.bed"
    bed = bedtools.GenomeCov(filename, f"{test_dir}/data/genbank/unicycler.gbk", chunksize=6000)
    chrom = bed.chr_list[0]
    chrom.run(501, k=2, circular=True)
    chrom.get_summary()
    chrom.get_rois()


def test_annotate_with_gff3():
    filename = f"{test_dir}/data/bed/subsample.bed"
    bed_gbk = bedtools.GenomeCov(filename, annotation_file=f"{test_dir}/data/genbank/subsample.gb")
    chrom = bed_gbk.chr_list[0]
    chrom.run(10001, k=2)
    chrom.get_summary()
    rois_gbk = chrom.get_rois()
    bed_gff = bedtools.GenomeCov(filename, annotation_file=f"{test_dir}/data/gff/subsample.gff3")
    chrom = bed_gff.chr_list[0]
    chrom.run(10001, k=2)
    chrom.get_summary()
    rois_gff = chrom.get_rois()
    # the only diff seems to be the index number
    assert rois_gff.df.equals(rois_gbk.df.reset_index(drop=True))
