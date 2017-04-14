# -*- coding: utf-8 -*-
"""
Created on Thu Jan 19 01:06:04 2017

@author: Shaurita
"""

from Bio.Phylo.PAML import codeml
import os
import csv

cml = codeml.Codeml()

# Set codeml options
cml.set_options(verbose = 0)
cml.set_options(CodonFreq = 2)
cml.set_options(cleandata = 1)
cml.set_options(fix_blength = 0)
cml.set_options(NSsites=[0, 1, 2])
cml.set_options(fix_omega = 0)
cml.set_options(clock = 1)
cml.set_options(ncatG = 2)
cml.set_options(runmode = 0)
cml.set_options(fix_kappa = 0)
cml.set_options(fix_alpha = 1)
cml.set_options(Small_Diff = 5e-7)
cml.set_options(method = 1)
cml.set_options(Malpha = 0)
cml.set_options(aaDist = 0)
cml.set_options(RateAncestor = 0)
cml.set_options(icode = 0)
cml.set_options(alpha = 0.0)
cml.set_options(seqtype = 1)
cml.set_options(omega = 0.4)
cml.set_options(getSE = 0)
cml.set_options(noisy = 3)
cml.set_options(Mgene = 0)
cml.set_options(kappa = 2)
cml.set_options(model = 1)
cml.set_options(ndata = 66)

##### PROGRAM #####
# Directory of phyml output files
home = '/work5/r2295/bin/Orthologs-Project/T1_Aligned'

# Read a list of gene names from a .csv file.
genes = open('tier1genes.csv')  # 1st column - List of genes used
file1 = csv.reader(genes)

# Let me know where I am right before I start the loop.
print("\n" + "The current working directory is " +
      os.getcwd() + (2 * "\n"))  # Print current working directory
Gene_count = 0


# For loop to iterate through gene list and create directories & run programs
for Gene in file1:  # Loop to create trees for each gene related aligned file
    Gene_count = Gene_count + 1

    # Create directories for PAML files
    gd = home + "/" + str(Gene[0])

    # Working directory
    wd = "./pamlout"

    os.chdir(gd)
    print(os.getcwd())

    # Run the codeml program
    cml.alignment = str(Gene[0]) + "_aligned.phy"
    cml.tree = str(Gene[0]) + "_aligned.phy_phyml_tree.txt"
    cml.out_file = str(Gene[0]) + "_tree.out"
    cml.run(verbose = True, command = "codeml")

    results = codeml.read(cml.out_file)