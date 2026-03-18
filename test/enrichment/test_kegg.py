import os
from io import BytesIO
from unittest.mock import MagicMock, patch

import pandas as pd
import pytest
from PIL import Image

from sequana.enrichment.kegg import KEGGPathwayEnrichment

from . import test_dir


# Minimal KGML for eco04122 with two gene entries
_KGML_STUB = """\
<?xml version="1.0"?>
<!DOCTYPE pathway SYSTEM "https://www.genome.jp/kegg/xml/KGML_v0.7.2_.dtd">
<pathway name="path:eco04122" org="eco" number="04122"
         title="Sulfur relay system"
         image="https://www.genome.jp/kegg/pathway/eco/eco04122.png"
         link="https://www.kegg.jp/kegg-bin/show_pathway?eco04122">
    <entry id="1" name="eco:b2530" type="gene"
           link="https://www.kegg.jp/dbget-bin/www_bget?eco:b2530">
        <graphics name="iscS" fgcolor="#000000" bgcolor="#BFFFBF"
                  type="rectangle" x="200" y="150" width="46" height="17"/>
    </entry>
    <entry id="2" name="eco:b3470" type="gene"
           link="https://www.kegg.jp/dbget-bin/www_bget?eco:b3470">
        <graphics name="tusA" fgcolor="#000000" bgcolor="#BFFFBF"
                  type="rectangle" x="300" y="200" width="46" height="17"/>
    </entry>
</pathway>
"""


def _make_png_bytes(width=800, height=600):
    """Return the bytes of a simple white PNG image."""
    buf = BytesIO()
    Image.new("RGB", (width, height), color="white").save(buf, format="PNG")
    buf.seek(0)
    return buf.getvalue()


def _mock_requests_get(url, **kwargs):
    response = MagicMock()
    response.raise_for_status = MagicMock()
    if url.endswith("/image"):
        response.content = _make_png_bytes()
    elif url.endswith("/kgml"):
        response.text = _KGML_STUB
    return response


def test_annotate_pathway_image(tmpdir):
    """_annotate_pathway_image should produce a PNG using KGML gene positions."""
    outpng = str(tmpdir.join("eco04122.png"))
    # Colors: b2530 in orange (up), b3470 in blue with white text (down)
    colors = {"b2530": "#FF6600", "b3470": "#0000CC,white"}

    with patch("sequana.enrichment.kegg.requests.get", side_effect=_mock_requests_get):
        ke = object.__new__(KEGGPathwayEnrichment)
        ke._annotate_pathway_image("eco04122", colors, outpng)

    assert os.path.exists(outpng)
    img = Image.open(outpng)
    assert img.format == "PNG"
    width, height = img.size
    assert width == 800 and height == 600

    # Verify that the gene boxes were painted with the expected colors.
    # b2530: center (200,150), width 46, height 17 → box (177,141)-(223,158)
    # Check top-left corner pixel of box (should not be covered by text)
    assert img.getpixel((178, 143))[:3] == (255, 102, 0), "b2530 box corner should be orange (#FF6600)"
    # b3470: center (300,200) → box (277,191)-(323,208)
    assert img.getpixel((278, 193))[:3] == (0, 0, 204), "b3470 box corner should be dark blue (#0000CC)"


@pytest.mark.xfail(reason="connection issue")
def test_ke(tmpdir):
    up = pd.read_csv(f"{test_dir}/data/ecoli_up_gene.csv")
    down = pd.read_csv(f"{test_dir}/data/ecoli_down_gene.csv")
    up = list(up.Name)
    down = list(down.Name)
    gene_lists = {"up": up, "down": down, "all": up + down}

    # used genes should be all genes used in the analyses. We do not have it in this example,
    # so we just add up and down
    ke = KEGGPathwayEnrichment(
        gene_lists, "eco", preload_directory=f"{test_dir}/data/kegg_pathways/", used_genes=up + down
    )

    # set background ourself
    ke = KEGGPathwayEnrichment(gene_lists, "eco", preload_directory=f"{test_dir}/data/kegg_pathways/", background=10000)

    # let use use the kegg genes (no background, no used_genes)
    ke = KEGGPathwayEnrichment(
        gene_lists,
        "eco",
        preload_directory=f"{test_dir}/data/kegg_pathways/",
    )

    with pytest.raises(ValueError):
        ke.barplot("dummy")

    ke.barplot("up")
    ke.barplot("down")
    ke.plot_genesets_hist()
    ke.scatterplot("down")
    assert ke.find_pathways_by_gene("moaA")
    ke.barplot_up_and_down()

    # save one pathway (just one to speed up things)
    outpng = tmpdir.join("test.png")
    df = pd.read_csv(f"{test_dir}/data/ecoli_all_gene.csv", index_col=0)
    ke.save_pathway("eco04122", df, filename=outpng)

    # save all pathways (same as input)
    path = tmpdir.mkdir("pathways_tmp")
    ke.save_pathways(str(path))

    ke.save_project("TEST", outdir=path)
