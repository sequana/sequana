from sequana.compare import RNADiffCompare

from . import test_dir

def test_rnadiff_volcano():

    PATH1 = f"{test_dir}/data/rnadiff/rnadiff_bowtie.csv"
    PATH2 = f"{test_dir}/data/rnadiff/rnadiff_salmon.csv"
    from sequana.rnadiff import RNADiffTable
    r1 = RNADiffTable(PATH1)
    r2 = RNADiffTable(PATH2)
    c = RNADiffCompare(r1, r2)


    c.plot_volcano()

    for mode in ['down', 'all', 'up']:
        c.plot_common_major_counts(mode)
        c.plot_jaccard_distance(mode)

    c.plot_foldchange()
    c.plot_volcano_differences() 
    c.plot_venn_all()
    c.plot_venn_up()
    c.plot_venn_down()
    #c.plot_corrplot_counts_raw()
    #c.plot_corrplot_counts_normed()
    #assert c.summary() == {"up1": 1295,   'up2': 56,
    #     'down1': 1325,
    #     'down2': 163,
    #     'common_down_r1_r2': 112,
    #     'common_up_r1_r2': 32}

