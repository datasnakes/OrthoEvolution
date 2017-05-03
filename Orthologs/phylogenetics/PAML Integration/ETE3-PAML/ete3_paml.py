# -*- coding: utf-8 -*-
"""
Date created: Wed Mar 15 15:14:05 2017
Author: S. Hutchins

Script description: Integrate ETE3 package with PAML & SLR on the MCSR.

"""
# Modules used
from ete3 import EvolTree
import os


# Import the newick tree
tree = EvolTree("ECP_EDN_15.nw")

# Import the alignment.
tree.link_to_alignment ("ECP_EDN_15.fasta")

# Set a working directory
tree.workdir = ""

# Set the binpath of the codeml binary
tree.execpath = r"C:\Users\shutchins2\Desktop\Software & Executables\paml4.9c\bin"

tree.run_model ('M1.example')
#x = tree.get_evol_model("M2")

## Run the codeml models
#models = ["M0", "M1", "M2"]
#for m in models:
#    tree.workdir = "PAML-Test"
#    tree.run_model = m
#    print("The " + m + " model ran successfully.")