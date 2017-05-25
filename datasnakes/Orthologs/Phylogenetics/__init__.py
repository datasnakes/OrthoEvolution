"""Phylogenetic tools part of the Orthologs Package"""

import warnings

from Bio import AlignIO
from datasnakes.Orthologs import OrthologsWarning

# Ignore the warning in this init script.
warnings.simplefilter('ignore', OrthologsWarning)

# Initialize the modules
from datasnakes.Orthologs.Phylogenetics import PhyML, ETE3PAML, Phylip, PamlTest
from datasnakes.Orthologs.Phylogenetics.PhyloTree import TreeViz

# Add a new module


class RelaxPhylip(object):
    """Convert the a multiple sequence alignment file to
    relaxed-phylip format.
    """
    def __init__(inputfile, outputfile):
        """Fasta to Relaxed Phylip format."""
        AlignIO.convert(inputfile, "fasta",
                        outputfile, "phylip-relaxed")


# Make this explicit, then they show up in the API docs
__all__ = ("ETE3PAML",
           "PamlTest",
           "PhyML",
           "TreeViz",
           "Phylip",
           "RelaxPhylip",
           )
