"""Phylogenetic Tools part of the Orthologs Package"""

import warnings
from Orthologs import OrthologsWarning

# Ignore the warning in this init script.
warnings.simplefilter('ignore', OrthologsWarning)

from Orthologs.Phylogenetics.PAML import ETE3PAML, PamlTest
from Orthologs.Phylogenetics.PhyML import PhyML
from Orthologs.Phylogenetics.Phylip import Phylip
from Orthologs.Phylogenetics.PhyloTree import TreeViz

# Make this explicit, then they show up in the API docs
__all__ = ("ETE3PAML",
           "PamlTest",
           "PhyML",
           "TreeViz",
           "Phylip",
)
