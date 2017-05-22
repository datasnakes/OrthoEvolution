# -*- coding: utf-8 -*-
"""
File Name:  TreeViz.py
Description: UNDER DEVELOPMENT!!!

Author: S. Hutchins
Date created: Wed Mar  1 14:41:34 2017
Project Name: Orthologs Project
"""
from Bio import Phylo

class TreeViz(object):
    def drawtree(path2tree, treeformat='newick'):
        """Import a newick formatted tree and visualize it."""
        tree = Phylo.read(path2tree, treeformat)
        Phylo.draw(tree)
