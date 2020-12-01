from sequana.find_motif import FindMotif
from sequana import sequana_data


def test_FindMotif():
    fm = FindMotif()
    fm.find_motif_bam(sequana_data("test_measles.bam"), 
        motif="CAG", local_threshold=5, global_threshold=10)

    fm.find_motif_from_sequence("CAG"*10, 
        motif="CAG", local_threshold=5)

    fm.find_motif_fasta(sequana_data('measles.fa'), 'CAG')

    fm.plot_alignment(sequana_data("test_measles.bam"), "CAG", window=30)

