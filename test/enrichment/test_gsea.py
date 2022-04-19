from sequana.enrichment.gsea import GSEA
import pytest
from . import test_dir



def test_gsea(tmpdir):
    # species is optional for now
    gs = GSEA("ecoli")
    gs.gene_sets = {'eco00592': ['pldA', 'fadI', 'fadA'],
        'eco04122': ['iscS', 'sseA', 'tusA', 'tusB', 'tusC', 'tusD', 'tusE', 'mnmA',
                     'moaE', 'moaA', 'moaD', 'ynjE', 'moeB', 'moaC', 'mog', 'moaB',
                     'thiI', 'thiS', 'thiF'],
        'eco00290': ['tdcB', 'ilvA', 'leuC', 'leuD', 'leuB', 'ilvI', 'ilvB', 'ilvH',
                     'ilvN', 'ilvM', 'ilvC', 'ilvD', 'ilvE', 'avtA', 'alaA', 'leuA']}
    gs.compute_enrichment(['tdcB', 'ilvA', 'leuC', 'leuD', 'leuB'], background=4000)
