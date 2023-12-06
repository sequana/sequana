from sequana import contigs

from . import test_dir


def test_all():
    filename = f"{test_dir}/data/fasta/test_contigs_ex1.fasta"
    c = contigs.Contigs(filename)
    c.df
    c.stats()
    c.hist_plot_contig_length()
    c.plot_contig_length_vs_GC()
    c.plot_contig_length_vs_nreads()

    c.plot_scatter_contig_length_vs_nreads_cov()
    c.plot_scatter_contig_length_vs_nreads_cov(min_nreads=0)
    c.plot_scatter_contig_length_vs_nreads_cov(min_nreads=0, logx=False)
    c.plot_scatter_contig_length_vs_nreads_cov(min_nreads=0, logy=False)
    c.scatter_length_cov_gc(logx=True, logy=True)


def test_spades():
    filename = f"{test_dir}/data/fasta/test_contigs_spades.fasta"
    c = contigs.Contigs(filename)
    c.df
    c.stats()
    c.hist_plot_contig_length()
    c.plot_contig_length_vs_GC()
    c.plot_contig_length_vs_nreads()
    c.plot_scatter_contig_length_vs_nreads_cov()
    c.plot_scatter_contig_length_vs_nreads_cov(min_nreads=0)
    c.plot_scatter_contig_length_vs_nreads_cov(min_nreads=0, logx=False)
    c.plot_scatter_contig_length_vs_nreads_cov(min_nreads=0, logy=False)
    c.scatter_length_cov_gc(logx=True, logy=True)
