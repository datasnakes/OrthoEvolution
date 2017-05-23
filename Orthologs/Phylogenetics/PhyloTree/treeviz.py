"""Import a newick formatted tree txt file and view it."""
from Bio import Phylo
from Orthologs import OrthologsDevelopmentWarning
import warnings

# Warn users about this module
warnings.warn('Orthologs.Phylogenetics.PhyloTree is still under development and '
              'may undergo significant changes prior to its official release.',
              OrthologsDevelopmentWarning)


class TreeViz(object):
    """Vizualize a newick tree."""

    def __init__(self, path2tree):
        """Import the path to the tree."""
        self.path2tree = path2tree

    def drawtree(self, treeformat='newick'):
        """Import a newick formatted tree and visualize it."""
        tree = Phylo.read(self.path2tree, treeformat)
        Phylo.draw(tree)
