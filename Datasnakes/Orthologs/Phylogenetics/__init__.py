"""Phylogenetic tools part of the Orthologs Package"""
import warnings
from Bio import AlignIO
from Datasnakes.Orthologs import OrthologsWarning

# Ignore the warning in this init script.
warnings.simplefilter('ignore', OrthologsWarning)

# Initialize the modules
from Datasnakes.Orthologs.Phylogenetics.PAML import ETE3PAML
from Datasnakes.Orthologs.Phylogenetics.PhyloTree import TreeViz
from Datasnakes.Orthologs.Phylogenetics import PhyML
from Datasnakes.Orthologs.Phylogenetics import Phylip
<<<<<<< HEAD
from Datasnakes.Orthologs.Phylogenetics.IQTree import IQTreeCommandline

# Ignore the warning in this init script.
warnings.simplefilter('ignore', OrthologsWarning)

=======
from Datasnakes.Orthologs.Phylogenetics import IQTree
>>>>>>> b4e6bb4ddfa7bc087f1fae5a8844594a6a6198c4

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
