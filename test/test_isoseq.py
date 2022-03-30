from sequana.isoseq import SIRV, SIRVReference
from easydev import md5, TempFile

from . import test_dir
def test_SIRVReference():

    with TempFile() as fh:
        data = f"{test_dir}/data/xls/test_sirv.xls"
        ss = SIRVReference()
        ss.from_excel(data)
        ss.to_fasta(fh.name)
        assert md5(fh.name) == "fa8a8b2d6cebf69bb84f5241592ac5f2"

def test_SIRV():

    with TempFile() as fh:
        data = f"{test_dir}/data/xls/test_sirv.xls"
        ss = SIRVReference()
        ss.from_excel(data)
        ss.to_fasta(fh.name)

        sirv = SIRV(fh.name)
        assert sirv.group_lengths == {'SIRV1': 7, 'SIRV2': 6, 'SIRV3': 11, 
            'SIRV4': 7, 'SIRV5': 12, 'SIRV6': 18, 'SIRV7': 7}
        assert sum(sirv.SIRV.lengths) == 75469


def test_mapped_sirv():
    from sequana import isoseq
    bam = f"{test_dir}/data/bam/test_isoseq_lq_sirv.bam"
    sirv = f"{test_dir}/data/fasta/SIRV.fa"
    mul = isoseq.PacbioIsoSeqMultipleIsoforms(sirv)
    mul.read(bam, "s1")
    mul.plot_bar()
    mul.plot_bar_grouped()
    mul.filter_raw_data()
    mul.plot_corr()
    list(mul.spikes_found().s1) == [1,1,1,1,1,1]
