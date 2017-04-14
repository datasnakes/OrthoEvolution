# -*- coding: utf-8 -*-
"""
Last updated on November 7, 2016

@author: Shaurita D. Hutchins

Edited and updated for use on the MCSR.
"""

# Integrate PhyML

# List of modules used
from Bio import Phylo
from Bio.Phylo.Applications import PhymlCommandline
import sys
from Bio import AlignIO

# Use the phyml executable file
phyml_exe = None

# This is mainly intended for windows use or use with an executable file
exe_name = "PhyML-3.1_win32.exe" if sys.platform == "win32" else "phyml"
phyml_exe = exe_name

# Convert the file to relaxed-phylip format
AlignIO.convert("egfr-family.aln", "fasta", "HTR1E_aligned.phy", "phylip-relaxed")

# Create the command & run phyml
# Input a phylip formatted alignment file and describe the datatype ('nt' or 'aa')
run_phyml = PhymlCommandline(phyml_exe, input='HTR1E_aligned.phy', datatype='nt')
print(run_phyml)
out_log, err_log = run_phyml()
