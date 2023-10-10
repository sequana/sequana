from sequana.rnadiff import RNADiffResults, RNADiffAnalysis, RNADesign, RNADiffTable
from . import test_dir
import pytest

def test_design():
    d = RNADesign(f"{test_dir}/data/rnadiff/design.csv")
    assert d.comparisons == [('Complemented_csrA', 'Mut_csrA'),  ('Complemented_csrA', 'WT'),  ('Mut_csrA', 'WT')]
    assert d.conditions ==  ['Complemented_csrA', 'Mut_csrA', 'WT']

    d = RNADesign(f"{test_dir}/data/rnadiff/design.csv", reference="WT")
    assert d.comparisons == [('Complemented_csrA', 'WT'),  ('Mut_csrA', 'WT')]
    assert d.conditions ==  ['Complemented_csrA', 'Mut_csrA', 'WT']


    try:
        d = RNADesign(f"{test_dir}/data/rnadiff/design_wrong.csv", reference="WT")
        assert False
    except (KeyError, SystemExit):
        assert True

    try:
        d = RNADesign(f"{test_dir}/data/rnadiff/design.csv", condition_col="TEST")
        assert False
    except SystemExit:
        assert True


@pytest.mark.xfail(reason="too slow or service may be down")
def test_rnadiff_onefolder():

    # Featurecounts are saved in sequana/resources/testing/rnadiff/rnadiff_onecond_ex1
    # generated from Featurecount of the file to be found in 
    # sequana/resources/testing/featurecounts/featurecounts_ex1

    counts = f"{test_dir}/data/rnadiff/rnadiff_onecond_ex1/counts.csv"
    design = f"{test_dir}/data/rnadiff/rnadiff_onecond_ex1/design.csv"
    gff = f"{test_dir}/data/rnadiff/rnadiff_onecond_ex1/Lepto.gff"

    # test minimum_mean_reads_per_condition_per_gene
    an = RNADiffAnalysis(counts, design, 
            condition="condition", comparisons=[("Complemented_csrA", "WT")], 
            fc_feature="gene", fc_attribute="ID", gff=gff, minimum_mean_reads_per_gene=1)

    # test wrong comparison
    try:
        an = RNADiffAnalysis(counts, design, 
            condition="condition", comparisons=[("Complemented_csrA", "W")], 
            fc_feature="gene", fc_attribute="ID", gff=gff, minimum_mean_reads_per_gene=1)
        assert False
    except SystemExit:
        assert True

    # test wrong ref
    try:
        an = RNADiffAnalysis(counts, design, 
            condition="condition", comparisons=[("Complemented_csrA", "WT"),], 
            fc_feature="gene", fc_attribute="ID", gff=gff, reference="WRONG")
        assert False
    except ValueError:
        assert True

    # test minimum_mean_reads_per_condition_per_gene
    an = RNADiffAnalysis(counts, design, 
            condition="condition", comparisons=[("Complemented_csrA", "WT")], 
            fc_feature="gene", fc_attribute="ID", gff=gff, minimum_mean_reads_per_condition_per_gene=1)
    an
    print(an)

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


def test_rnadiffTable():
    table = f"{test_dir}/data/rnadiff/rnadiff_onecond_1/tables/B3789-v1.surexpvsref.complete2.xls"

    # 
    r = RNADiffTable(table, sep="\t", shrinkage=False)

    r = RNADiffTable(table, sep="\t", shrinkage=True)
    r.alpha = 0.04
    assert r.alpha==0.04
    r.log2_fc = 1.1
    assert r.log2_fc ==1.1
    r.summary()
    r.plot_volcano(plotly=True)

    # FIXME: fails on py3.9 and 3.10. success on 3.8
    # r.plot_volcano(add_broken_axes=True)
    r.plot_volcano(add_broken_axes=False)
    r.plot_pvalue_hist()
    r.plot_padj_hist()

def test_rnadiffResults():    

    rnadiff = f"{test_dir}/data/rnadiff/rnadiff_0.15.4"
    r = RNADiffResults(rnadiff)
    r.alpha = 0.04
    r.log2_fc = 1
    #r.plot_dispersion()
    #comp = list(r.comparisons.keys())[0]
    #r.heatmap_vst_centered_data(comp)
    r.plot_count_per_sample()
    r.plot_feature_most_present()
    r.plot_most_expressed_features()
    r.plot_percentage_null_read_counts()
    r.plot_pca(plotly=True, n_components=3)
    r.plot_pca(plotly=False)
    r.plot_mds()
    r.plot_upset()
