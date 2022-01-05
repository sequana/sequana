from sequana.pacbio import BAMSimul, PBSim, PacbioSubreads, Barcoding, CCS
from sequana import sequana_data
from easydev import TempFile


def test_pacbio():
    b = PacbioSubreads(sequana_data("test_pacbio_subreads.bam"))
    assert len(b) == 130
    b.df
    #assert b.nb_pass[1] == 130


    with TempFile() as fh:
        b.filter_length(fh.name, threshold_min=500)

    print(b)   #  check length

    assert b.stats['mean_GC'] > 62.46
    assert b.stats['mean_GC'] < 65.47

    b.summary()

    b = PacbioSubreads(sequana_data("test_pacbio_subreads.bam"))

    # test hist_snr from scratch
    b._df = None
    b.hist_snr()

    # test hist_len from scratch
    b._df = None
    b.hist_read_length()
    b.hist_nb_passes()
    b.get_mean_nb_passes()

    # test from scratch
    b._df = None
    b.hist_GC()

    # test from scratch
    b._df = None
    b.plot_GC_read_len()

    # test from scratch
    b._df = None

    with TempFile() as fh:
        b.to_fasta(fh.name, threads=1)
    with TempFile() as fh:
        b.to_fastq(fh.name, threads=1)
    with TempFile() as fh:
        b.save_summary(fh.name)


def test_pacbio_stride():
    b = PacbioSubreads(sequana_data("test_pacbio_subreads.bam"))
    with TempFile() as fh:
        b.stride(fh.name, stride=2)
    with TempFile() as fh:
        b.stride(fh.name, stride=2, random=True)


def test_pacbio_random():
    b = PacbioSubreads(sequana_data("test_pacbio_subreads.bam"))

    with TempFile() as fh:
        b.random_selection(fh.name, nreads=10)

    with TempFile() as fh:
        b.random_selection(fh.name,  expected_coverage=10,
            reference_length=10000)


def test_bamsim():
    filename = sequana_data("test_pacbio_subreads.bam")
    b = BAMSimul(filename)
    b.df
    b.hist_read_length()
    b.hist_GC()
    b.plot_GC_read_len()
    with TempFile() as fh:
        b.filter_length(fh.name, threshold_min=500)
    with TempFile() as fh:
        mask = [True for this in range(len(b))]
        b.filter_bool(fh.name, mask)


def test_pbsim():
    filename = sequana_data("test_pacbio_subreads.bam")
    ss = PBSim(filename, filename)
    with TempFile() as fh:
        ss.run(bins=100, step=50, output_filename=fh.name)
    from pylab import close
    close()


def test_barcoding():
    data = sequana_data("test_pacbio_barcode_report.csv")
    bc = Barcoding(data)

    import tempfile
    with tempfile.TemporaryDirectory() as tempdir:
        bc.plot_and_save_all(directory=tempdir)


