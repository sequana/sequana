from sequana.compare import RNADiffCompare
from sequana import sequana_data


def test_rnadiff_volcano():


    PATH1 = sequana_data("rnadiff/rnadiff_twocond_ex1")
    PATH2 = sequana_data("rnadiff/rnadiff_twocond_ex2")
    gff = sequana_data("rnadiff/rnadiff_onecond_ex1/Lepto.gff")

    from sequana.rnadiff import RNADiffResults
    r1 = RNADiffResults(PATH1, PATH1 + "/../design_twocond.csv",  
        fc_feature="gene", fc_attribute="ID", gff=gff)
    r2 = RNADiffResults(PATH2, PATH1 + "/../design_twocond.csv",  
        fc_feature="gene", fc_attribute="ID", gff=gff)
    c = RNADiffCompare(r1, r2)

    c.comparison = 'Complemented_csrA_vs_Mut_csrA'

    c.plot_volcano()

    for mode in ['down', 'all', 'up']:
        c.plot_common_major_counts(mode)
        c.plot_jaccard_distance(mode)
       
    c.plot_foldchange()
    c.plot_volcano_differences() 
    c.plot_venn_all()
    c.plot_venn_up()
    c.plot_venn_down()
    c.plot_corrplot_counts_raw()
    c.plot_corrplot_counts_normed()
    #assert c.summary() == {"up1": 1295,   'up2': 56,
    #     'down1': 1325,
    #     'down2': 163,
    #     'common_down_r1_r2': 112,
    #     'common_up_r1_r2': 32}

