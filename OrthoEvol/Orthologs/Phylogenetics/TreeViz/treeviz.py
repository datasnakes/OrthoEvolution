"""Import a newick formatted tree txt file and view it."""
import warnings

from Bio import Phylo
from ete3 import Tree

from OrthoEvol.Orthologs import OrthologsDevelopmentWarning


class TreeViz(object):
    """Tools that allow visualization of a newick formatted tree."""

    def __init__(self, path, tree_format='newick'):
        """Initialize the class.

        :param path:  The path to your tree file.
        :type path: str
        :param tree_format:  The format of the tree, default value = 'newick'
        :type path: tree_format
        """
        # Warn users about this module
        warnings.warn('This module is still under development and '
                      'may undergo significant changes prior to its official '
                      'release.', OrthologsDevelopmentWarning)
        self.path = path
        self.tree_format = tree_format
        # Read the tree
        self.tree = self.read_tree(path=path, tree_format=tree_format)

    def read_tree(self, path, tree_format):
        """Read the phylogenetic tree.

        :param path: The path to your tree file.
        :type path: str
        :param tree_format: The format of the tree, default value = 'newick'
        :type tree_format: str
        """
        tree = Phylo.read(file=self.path, format=self.tree_format)
        return tree

    def draw_tree(self, drawing_type="default"):
        """Import a newick formatted tree and visualize it."""
        if drawing_type == "ascii":
            Phylo.draw_ascii(self.tree)
        elif drawing_type == "graphviz":
            Phylo.draw_graphviz(self.tree)
        elif drawing_type == "default":
            Phylo.draw(self.tree)
