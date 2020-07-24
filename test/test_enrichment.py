from sequana.enrichment import KeggPathwayEnrichment
from sequana import sequana_data
from easydev import TempFile



def test_ke():
    RNADIFF_DIR = sequana_data("rnadiff") + "/rnadiff_onecond_1"
    ke = KeggPathwayEnrichment(RNADIFF_DIR, "eco", log2_fc=1)
    ke.barplot(ke.enrichment['down'])
    ke.barplot(ke.enrichment['up'])
    ke.barplot(ke.enrichment['down'])
    ke.plot_genesets_hist()
    ke.scatterplot(ke.enrichment['down'])
    assert ke.find_pathways_by_gene("moaA")
    assert ke.find_pathways_by_gene("moaA", match="exact")

    with TempFile(suffix=".png") as fout:
        ke.save_pathway(ke.kegg.pathwayIds[0].replace(":", "").replace("path",
""), filename=fout.name)


def test_panther():
    from sequana.enrichment import PantherEnrichment
    RNADIFF_DIR = sequana_data("rnadiff") + "/rnadiff_onecond_1"
    pe = PantherEnrichment(RNADIFF_DIR, log2_fc_threshold=0)
    pe = PantherEnrichment(RNADIFF_DIR, fc_threshold=1)

    pe.compute_enrichment(pe.mygenes_up,  83333, ontologies=["GO:0003674"],
        correction="bonferroni")


    # too slow or fails ith uniprot 404 from time to time
    #pe.get_functional_classification(pe.mygenes_up, 83333)


    ontologies = ["GO:0003674", "GO:0008150", "GO:0005575"]



    pe.compute_enrichment(pe.mygenes_up,  83333, ontologies=ontologies)

    df = pe.plot_go_terms(ontologies=ontologies, compute_levels=False)
    df = pe.plot_go_terms(ontologies=ontologies, compute_levels=False, log=True)
    df = pe.plot_go_terms(ontologies=ontologies, compute_levels=False, log=True,
        include_negative_enrichment=True)

    pe.plot_piechart(df) 

    #with TempFile(suffix=".png") as fout:
    #    pe.save_chart(df, filename=fout.name)
    # too slow
