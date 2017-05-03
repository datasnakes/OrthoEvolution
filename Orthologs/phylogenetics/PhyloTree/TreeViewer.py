# -*- coding: utf-8 -*-
"""
Date created: Wed Mar  1 14:41:34 2017
Author: S. Hutchins

Script description:

"""
from Bio import Phylo

tree = Phylo.read('apol1_tree.txt', 'newick')
Phylo.draw(tree)
