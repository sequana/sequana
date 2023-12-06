import os
from collections import defaultdict

from sequana.fasta_gff_correction import FastaGFFCorrection

from . import test_dir


def test_fasta_gff_correction(tmpdir):
    # test with a custom fasta
    in_fasta = f"{test_dir}/data/fasta_gff_correction/test.fa"
    in_gff = f"{test_dir}/data/fasta_gff_correction/test.gff3"
    in_vcf = f"{test_dir}/data/fasta_gff_correction/test.vcf"

    out_fasta = tmpdir.join("test.fa")
    out_gff = tmpdir.join("test.gff3")

    f = FastaGFFCorrection(in_fasta, in_vcf)
    f.fix_and_save_fasta(out_fasta)
    f.fix_and_save_gff(in_gff, out_fasta, out_gff)

    # we can read the output gff and output sequence back and count the stop/start codons

    assert f.get_all_start_codons(out_fasta, out_gff, "AE000666_1", strand="+") == defaultdict(
        int, {"ATG": 50, "TTG": 17, "GTG": 13}
    )

    assert f.get_all_stop_codons(out_fasta, out_gff, "AE000666_1", strand="+") == defaultdict(
        int, {"TGA": 29, "TAA": 31, "TAG": 20}
    )

    assert f.get_all_start_codons(out_fasta, out_gff, "AE000666_1", strand="-") == defaultdict(
        int, {"CAT": 43, "CAC": 9, "CAA": 12}
    )

    assert f.get_all_stop_codons(out_fasta, out_gff, "AE000666_1", strand="-") == defaultdict(
        int, {"TTA": 35, "TCA": 17, "CTA": 12}
    )
