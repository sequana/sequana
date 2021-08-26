from sequana.enrichment.mart import Mart
import pytest


@pytest.mark.xfail(reason="too slow or service may be down")
def test_mart():
    conv = Mart(dataset="mmusculus_gene_ensembl")
    # you could choose hsapiens_gene_ensembl for instance
    df = conv.query()
    df.set_index("ensembl_gene_id")
    #    conv.save(df)

