# -*- coding: utf-8 -*-
"""
Date created: Mon Mar 13 16:12:43 2017
Author: S. Hutchins

Script description: Quick phyml to paml test on one multiple sequence alignment.

"""
from Bio.Phylo.PAML import codeml
from Bio import AlignIO
from Bio.Phylo.Applications import PhymlCommandline
import sys

#------------------------------------------------------------------------------
# Convert the file to relaxed-phylip format
AlignIO.convert("HTR1A_aligned_cds_nucl.fasta", "fasta", "HTR1A_aligned.phy", "phylip-relaxed")

#------------------------------------------------------------------------------
# Run phyml to generate tree results

# Use the phyml executable file
phyml_exe = None

# This is mainly intended for windows use or use with an executable file
exe_name = "PhyML-3.1_win32.exe" if sys.platform == "win32" else "phyml"
phyml_exe = exe_name

# Create the command & run phyml
# Input a phylip formatted alignment file and describe the datatype ('nt' or 'aa')
run_phyml = PhymlCommandline(phyml_exe, input='HTR1A_aligned.phy', datatype='nt')
print(run_phyml)
out_log, err_log = run_phyml()

#------------------------------------------------------------------------------
# Run PAML. More specifically use the codeml program

cml = codeml.Codeml()


# Set PAML control file options
cml.set_options(verbose = 1)
cml.set_options(CodonFreq = 2)
cml.set_options(cleandata = 1)
cml.set_options(fix_blength = 0)
cml.set_options(NSsites=[0])
cml.set_options(fix_omega = 0)
cml.set_options(clock = 0)
cml.set_options(ncatG = 2)
cml.set_options(runmode = 0)
cml.set_options(fix_kappa = 0)
cml.set_options(fix_alpha = 1)
cml.set_options(Small_Diff = 5e-7)
cml.set_options(method = 0)
cml.set_options(Malpha = 0)
cml.set_options(aaDist = 0)
cml.set_options(RateAncestor = 0)
cml.set_options(icode = 0)
cml.set_options(alpha = 0.0)
cml.set_options(seqtype = 1)
cml.set_options(omega = 0.4)
cml.set_options(getSE = 0)
cml.set_options(noisy = 9)
cml.set_options(Mgene = 0)
cml.set_options(kappa = 2)
cml.set_options(model = 1)
cml.set_options(ndata = 1)

# Run PAML (codeml in this example)
# Directory of phyml input and output files
home = '/work5/r2295/bin/Orthologs-Project/Tests'

cml.read_ctl_file("codeml.ctl")
cml.tree = "HTR1A_aligned.phy_phyml_tree.txt"
cml.alignment = "HTR1A_aligned.phy"
cml.out_file = "output.txt"
cml.working_dir = home

cml.run(verbose=True, command = "codeml")
results = cml.read(cml.out_file)