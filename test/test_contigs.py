from sequana import contigs, sequana_data


def test_all():
    filename = sequana_data('test_contigs_ex1.fasta')
    c = contigs.Contigs(filename, filename)
    c.stats()
    c.plot_contig_length_vs_GC()
    #c.plot_contig_length_vs_nreads()

    c.hist_plot_contig_length()
    c.get_df()
