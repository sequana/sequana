from sequana.rnadiff import RNADiffResults, RNADiffAnalysis, RNADesign
from sequana import sequana_data


def test_design():
    d = RNADesign(sequana_data("rnadiff/design.csv"))
    assert d.comparisons == [('Complemented_csrA', 'Mut_csrA'),  ('Complemented_csrA', 'WT'),  ('Mut_csrA', 'WT')]
    assert d.conditions ==  ['Complemented_csrA', 'Mut_csrA', 'WT']

    d = RNADesign(sequana_data("rnadiff/design.csv"), reference="WT")
    assert d.comparisons == [('Complemented_csrA', 'WT'),  ('Mut_csrA', 'WT')]
    assert d.conditions ==  ['Complemented_csrA', 'Mut_csrA', 'WT']

def test_rnadiff_onefolder():



    # Featurecounts are saved in sequana/resources/testing/rnadiff/rnadiff_onecond_ex1

    # generated from Featurecount of the file to be found in 
    # sequana/resources/testing/featurecounts/featurecounts_ex1

    counts = sequana_data("rnadiff/rnadiff_onecond_ex1/counts.csv")
    design = sequana_data("rnadiff/rnadiff_onecond_ex1/design.csv")
    gff = sequana_data("rnadiff/rnadiff_onecond_ex1/Lepto.gff")

    an = RNADiffAnalysis(counts, design, 
            condition="condition", comparisons=[("Complemented_csrA", "WT")], 
            fc_feature="gene", fc_attribute="ID", gff=gff)
    an

    r = an.run()

    r.plot_count_per_sample()
    r.plot_percentage_null_read_counts()
    #r.plot_volcano()
    r.plot_pca()
    r.plot_mds()
    r.plot_isomap()
    r.plot_density()
    r.plot_boxplot_normeddata()
    r.plot_boxplot_rawdata()
    r.plot_dendogram()
    r.plot_dispersion()
    r.plot_feature_most_present()

    r.comparisons['Complemented_csrA_vs_WT'].plot_volcano()
    r.comparisons['Complemented_csrA_vs_WT'].plot_padj_hist()
    r.comparisons['Complemented_csrA_vs_WT'].plot_pvalue_hist()

    r.summary()
    r.alpha = 1
    r.log2_fc = 1
