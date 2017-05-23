"""
Import a newick formatted tree txt file and view it.
"""
from Bio import Phylo
from Orthologs import OrthologsDevelopmentWarning
import warnings

# Warn users about this module
warnings.warn('Orthologs.Phylogenetics.PhyloTree is still under development and '
              'may undergo significant changes prior to its official release.',
              OrthologsDevelopmentWarning)


class TreeViz(object):
    def __init__(self, treepath):
        path2tree = treepath
        return path2tree

    def drawtree(path2tree, treeformat='newick'):
        """Import a newick formatted tree and visualize it."""
        tree = Phylo.read(path2tree, treeformat)
        Phylo.draw(tree)
