import pytest
from genomeuploader.taxon_finder import TaxonFinder

def test_taxon_finder_bacteria():
    # unknown species
    finder = TaxonFinder('d__Bacteria;p__Pseudomonadota;c__Alphaproteobacteria;o__;f__;g__;s__')
    assert finder.taxid == 1913988
    assert finder.scientific_name == 'Alphaproteobacteria bacterium'

    # full lineage
    finder = TaxonFinder('d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Phocaeicola;s__Phocaeicola plebeius')
    assert finder.taxid == 310297
    assert finder.scientific_name == 'Phocaeicola plebeius'

    # species with sp.
    finder = TaxonFinder('d__Bacteria;p__Bacillota;c__Clostridia;o__Eubacteriales;f__Oscillospiraceae;g__Ruminococcus;s__')
    assert finder.taxid == 41978
    assert finder.scientific_name == 'Ruminococcus sp.'
