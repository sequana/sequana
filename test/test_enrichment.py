from sequana.enrichment import KeggPathwayEnrichment, PantherEnrichment
from sequana import sequana_data
from easydev import TempFile
import pandas as pd



def test_ke():
    up = pd.read_csv(sequana_data("enrichment/ecoli_up_gene.csv"))
    down = pd.read_csv(sequana_data("enrichment/ecoli_down_gene.csv")    )
    up = list(up.Name)
    down = list(down.Name)
    gene_lists = {'up': up, 'down': down, 'all': up +down}

    ke = KeggPathwayEnrichment(gene_lists, "lbi", log2_fc=0)
    ke.barplot('down')
    ke.barplot('up')
    ke.barplot('down')
    ke.plot_genesets_hist()
    ke.scatterplot('down')
    assert ke.find_pathways_by_gene("moaA")
    assert ke.find_pathways_by_gene("moaA", match="exact")

    with TempFile(suffix=".png") as fout:
        ke.save_pathway(ke.kegg.pathwayIds[0].replace(":", "").replace("path", ""), filename=fout.name)


def test_panther():
    up = pd.read_csv(sequana_data("enrichment/ecoli_up_gene.csv"))
    down = pd.read_csv(sequana_data("enrichment/ecoli_down_gene.csv")    )
    up = list(up.Name)
    down = list(down.Name)
    gene_lists = {'up': up, 'down': down, 'all': up +down}

    pe = PantherEnrichment(gene_lists, log2_fc_threshold=0, taxon=83333)
    pe = PantherEnrichment(gene_lists, fc_threshold=1, taxon=83333)

    pe.compute_enrichment(ontologies=["GO:0003674"],
        correction="bonferroni")


    # too slow or fails ith uniprot 404 from time to time
    #pe.get_functional_classification(pe.mygenes_up, 83333)


    ontologies = ["GO:0003674", "GO:0008150", "GO:0005575"]



    pe.compute_enrichment(  ontologies=ontologies)

    df = pe.plot_go_terms("up", ontologies=ontologies, compute_levels=False)
    df = pe.plot_go_terms("up", ontologies=ontologies, compute_levels=False, log=True)
    df = pe.plot_go_terms("up", ontologies=ontologies, compute_levels=False, log=True,
        include_negative_enrichment=True)

    pe.plot_piechart(df) 

    #with TempFile(suffix=".png") as fout:
    #    pe.save_chart(df, filename=fout.name)
    # too slow
