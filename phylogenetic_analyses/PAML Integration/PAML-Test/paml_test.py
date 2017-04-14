# -*- coding: utf-8 -*-
"""
Created on Thu Jan 19 12:55:18 2017

@author: shutchins2
"""

from Bio.Phylo.PAML import codeml

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
cml.tree = "ADRA1A_aligned.phy_phyml_tree.txt"
cml.alignment = "ADRA1A_aligned.phy"
cml.out_file = "output.txt"
cml.working_dir = home

cml.run(verbose = True, command = "/work5/r2295/bin/PAML/paml48/codeml")
results = codeml.read(cml.out_file)