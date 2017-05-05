# -*- coding: utf-8 -*-
"""
Date created: Wed Mar 15 15:14:05 2017
Author: S. Hutchins

Script description: Integration of ETE3 for using PAML's codeml. For this project,
we start my using the M1 model so it is the default.

"""
# Modules used
from ete3 import EvolTree


def ete3paml(gene, paml_path, workdir='data/paml-output/', model='M1'):
    # Import the newick tree
    tree = EvolTree("data/phyml-output/" + gene + "_PhyML/" + gene + "_tree.nw")

    # Import the alignment
    tree.link_to_alignment("data/clustal-output/" + gene + "_Aligned/" + gene + "_aligned_cds_nucl.fasta")

    tree.workdir = workdir

    # Set the binpath of the codeml binary
    tree.execpath = paml_path

    # Run the codeml model
    tree.run_model(model + '.' + gene)