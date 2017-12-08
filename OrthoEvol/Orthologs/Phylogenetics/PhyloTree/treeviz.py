"""Import a newick formatted tree txt file and view it."""
import warnings

from Bio import Phylo

from OrthoEvol.Orthologs import OrthologsDevelopmentWarning

# Warn users about this module
warnings.warn('This module is still under development and '
              'may undergo significant changes prior to its official release.',
              OrthologsDevelopmentWarning)


class TreeViz(object):
    """Tools that allow visualization of a newick formatted tree."""

    def __init__(self, path2tree, treeformat='newick'):
        """Import the path to the tree.

        :param path2tree:  Path to your tree file.
        :param treeformat:  (Default value = 'newick')
        """
        self.path2tree = path2tree
        self.treeformat = treeformat
        self.tree = Phylo.read(self.path2tree, self.treeformat)

    def drawtree(self):
        """Import a newick formatted tree and visualize it."""

        Phylo.draw(self.tree)
