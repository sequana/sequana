from sequana.sequence import DNA, RNA, Repeats

from . import test_dir

datafile = f"{test_dir}/data/fasta/measles.fa"


def test_dna():
    seq = "ACGTTTT"
    dna = DNA(seq)
    assert dna.sequence == seq
    assert dna.get_complement() == "TGCAAAA"
    assert dna.get_reverse() == "TTTTGCA"
    assert dna.get_reverse_complement() == "AAAACGT"
    dna.check()
    dna.stats()

    # inplace functions
    dna.reverse()
    dna.complement()
    dna.reverse_complement()

    dna = DNA("jjjj")
    try:
        dna.check()
        assert False
    except:
        assert True


    # read a file and tests the __len__ method
    dna = DNA(datafile)
    assert len(dna) == 15894
    dna.gc_content()
    # tests for the ORF and CDS functions
    assert dna.ORF_pos.shape == (1525, 6)
    dna.type_filter = "CDS"
    assert dna.ORF_pos.shape == (328, 6)
    dna.threshold = 20
    assert dna.ORF_pos.shape == (231, 6)
    dna.type_filter = "ORF"
    assert dna.ORF_pos.shape == (964, 6)
    dna.threshold = 0
    assert dna.ORF_pos.shape == (1525, 6)
    dna.hist_ORF_CDS_linearscale()
    dna.barplot_count_ORF_CDS_by_frame()

    # test occurences
    dna._data = "ACGTGGGGGTT"
    assert dna.get_occurences("GGG", False) == [4]
    assert dna.get_occurences("GGG", True) == [4, 5, 6]


def test_rna():
    rna = RNA("ACUG")


def test_repeats():
    rep = Repeats(datafile)
    rep.threshold = 11
    assert len(rep.begin_end_repeat_position) == 158
    rep.do_merge = True
    assert len(rep.begin_end_repeat_position) == 156
    rep.do_merge = False
    assert len(rep.begin_end_repeat_position) == 158
    assert rep.df_shustring.shape == (15888,2)
    assert rep.longest_shustring == 15

    # test histogram
    rep = Repeats(datafile)
    rep.threshold = 11
    rep.hist_length_repeats()



def test_gc_skew():

    dna = DNA(datafile)
    dna.window = 100
    dna.plot_all_skews()
