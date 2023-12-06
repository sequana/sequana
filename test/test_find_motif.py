from sequana.find_motif import FindMotif

from . import test_dir


def test_FindMotif():
    fm = FindMotif()
    fm.find_motif_bam(f"{test_dir}/data/bam/test_measles.bam", motif="CAG", local_threshold=5, global_threshold=10)

    fm.find_motif_from_sequence("CAG" * 10, motif="CAG", local_threshold=5)

    fm.find_motif_fasta(f"{test_dir}/data/fasta/measles.fa", "CAG")

    fm.plot_alignment(f"{test_dir}/data/bam/test_measles.bam", "CAG", window=30)
