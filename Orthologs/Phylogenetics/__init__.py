"""Phylogenetic Tools part of the Orthologs Package"""

import warnings
from Orthologs import OrthologsWarning

# Ignore the warning in this init script.
warnings.simplefilter('ignore', OrthologsWarning)

#from PAML.__init__ import ETE3PAML, PamlTest
#from PhyML.__init__ import PhyML, PhymlTest
#from Phylip.__init__ import Phylip
#from PhyloTree.__init__ import TreeViz
#
#
#
## Make this explicit, then they show up in the API docs
#__all__ = ("ETE3PAML",
#           "PamlTest",
#           "PhyML",
#           "PhymlTest",
#           "TreeViz",
#           "Phylip",
#)
