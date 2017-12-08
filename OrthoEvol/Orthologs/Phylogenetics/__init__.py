"""Phylogenetic tools part of the Orthologs Package."""
import warnings
from Bio import AlignIO
from OrthoEvol.Orthologs import OrthologsWarning

# Ignore the warning in this init script.
warnings.simplefilter('ignore', OrthologsWarning)

# Initialize the modules
from OrthoEvol.Orthologs.Phylogenetics.PAML import ETE3PAML
from OrthoEvol.Orthologs.Phylogenetics.PhyloTree import TreeViz
from OrthoEvol.Orthologs.Phylogenetics import PhyML
from OrthoEvol.Orthologs.Phylogenetics import Phylip
from OrthoEvol.Orthologs.Phylogenetics.IQTree import IQTreeCommandline

# Ignore the warning in this init script.
warnings.simplefilter('ignore', OrthologsWarning)


class RelaxPhylip(object):
    """Convert a multiple sequence alignment file to relaxed-phylip format."""

    def __init__(inputfile, outputfile):
        """Fasta to Relaxed Phylip format."""
        AlignIO.convert(inputfile, "fasta",
                        outputfile, "phylip-relaxed")


# Make this explicit, then they show up in the API docs
__all__ = ("ETE3PAML",
           "PhyML",
           "TreeViz",
           "RelaxPhylip",
           "Phylip",
           "IQTreeCommandline",
           )
