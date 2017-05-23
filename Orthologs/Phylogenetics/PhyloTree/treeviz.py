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

    def __init__(self):
        """Import the path to the tree."""
        self.path2tree

    def drawtree(self, treepath, treeformat='newick'):
        """Import a newick formatted tree and visualize it."""
        self.path2tree = treepath
        tree = Phylo.read(treepath, treeformat)
        Phylo.draw(tree)
