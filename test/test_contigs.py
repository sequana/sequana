from sequana import contigs
from . import test_dir

def test_all():
    filename = f"{test_dir}/data/fasta/test_contigs_ex1.fasta"
    c = contigs.Contigs(filename, filename)
    c.stats()
    c.plot_contig_length_vs_GC()
    #c.plot_contig_length_vs_nreads()

    c.hist_plot_contig_length()
    c.get_df()



def test_spades():
    filename = f"{test_dir}/data/fasta/test_contigs_spades.fasta"
    c = contigs.ContigsSpades(filename)
    c.hist_contig_length()
    c.plot_contig_length_vs_GC()
    c.scatter_length_cov_gc()
