# -*- coding: utf-8 -*-
"""
Created on Mon Jan  9 13:15:22 2017
@author: S. Hutchins

Interfacing with PAML via Biopython -> http://biopython.org/wiki/PAML

PAML (for Phylogenetic Analysis by Maximum Likelihood) is a package of
programs for phylogenetic analyses of DNA and protein sequences using maximum
likelihood.
"""

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
# Interfacing with PAML
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------


# List of modules used
from Bio.Phylo.PAML import codeml, baseml, yn00
from Bio.Phylo.PAML.chi2 import cdf_chi2
import csv
import pandas as pd
import os

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

    # Use the command line to call/run the executable programs
    #codeml = 'bin/codeml'
    #baseml = 'bin/baseml'
    #pamp = 'bin/pamp'

    #------------------------------------------------------------------------------
    # The Codeml Object
    #------------------------------------------------------------------------------
    """
    cml.print_options() - The codeml runtime options are stored in a dictionary
    object that is keyed by the option names. Options may be set by the
    set_option() function and their values may be retrieved by the get_option()
    function.

    cml.set_options(clock=1)
    cml.set_options(NSsites=[0, 1, 2])
    cml.set_options(aaRatefile="wag.dat")
    cml.get_option("NSsites")


    The control file to read is provided as an argument to the read_ctl_file()
    method, while the write_ctl_file() method writes to the Codeml object’s
    ctl_file attribute.

    cml = codeml.codeml()
    cml.read_ctl_file("codeml.ctl")
    cml.print_options()
    cml.ctl_file = "control2.ctl"
    cml.write_ctl_file()
    cml.run(verbose=True, parse=False, ctl_file=path_to_control_file, command=path_to_codeml_executable)
    results = codeml.read()


    Example for omega values for each branch:
    from Bio.Phylo.PAML import codeml
    results = codeml.read(paml_outputfile)
    print(results["NSsites"][0]["parameters"]["omega"])
    This gives you a list of omega (dn/ds) for each branch


    # The output below can be done individually
    cml = codeml.Codeml()
    cml.alignment = "align.phylip"
    cml.tree = "species.tree"
    cml.out_file = "results.out"
    cml.working_dir = "./scratch"

    When I get the amino acid sequences, I'll also use codeml with them.

    """

    # Create variables for the codeml commands
    cml = codeml.Codeml()

    # Read the example control file
    cml.read_ctl_file("/work5/r2295/bin/PAML/paml48/codonml.ctl")

    # Create an iniial control file
    cml.ctl_file = str(Gene[0]) + "_initialfile.ctl"
    cml.write_ctl_file()

    # Create an initial tree file
    cml.set_options(NSsites=[0])
    cml.set_options(aaRatefile=0)
    cml.set_options(CodonFreq=2)
    cml.set_options(model=0)
    cml.set_options(seqtype=1)
    cml.set_options(runmode=0)  # 0 = user tree,
    cml.set_options(ndata=66)
    cml.write_ctl_file()

    # Run the codeml program
    cml.alignment = str(Gene[0]) + "_aligned.phy"
    cml.tree = str(Gene[0]) + "_aligned.phy_phyml_tree.txt"
    cml.out_file = str(Gene[0]) + "_initialtree.out"
    cml.run(ctl_file=None, verbose=True, command="codeml")


    # Create a working control file
    cml.ctl_file = str(Gene[0]) + "_controlfile.ctl"
    cml.write_ctl_file()

    # Modify the parameters of the control file
    # This is for a free ratios test
    cml.set_options(NSsites=[0, 1 , 2])
    cml.set_options(model=1)
    cml.set_options(seqtype=1)
    cml.set_options(runmode=0)  # -2 is for pairwise
    cml.set_options(CodonFreq=2)
    cml.set_options(verbose=1)
    cml.write_ctl_file()

    # Run the codeml program
    cml.alignment = str(Gene[0]) + "_aligned.phy"
    cml.tree = str(Gene[0]) + "_initialtree.out"
    cml.out_file = str(Gene[0]) + "_results.out"
    cml.run(ctl_file=str(Gene[0]) + "_controlfile.ctl", verbose=True, command="codeml")


    #------------------------------------------------------------------------------
    # Baseml and Yn00
    #------------------------------------------------------------------------------
    """
    Baseml and Yn00 share the same methods and attributes as Codeml, and are thus
    used in the same manner. It should be noted, however, that Yn00 does not have
    a tree attribute, as yn00 does not require a tree file.

    It's best to use codeml for coding sequences, but I can also use baseml for
    the cds sequnces. I likely will do this at some point.

    """
#    bml = baseml.Baseml()
#    bml.print_options()
#
#    yn = yn00.Yn00()
#    yn.print_options()

    #------------------------------------------------------------------------------
    # chi2
    #------------------------------------------------------------------------------
    """
    The chi2 module offers an easy method to retrieve p-values from a Chi-squared
    cumulative distribution function for likelihood ratio tests, which are
    performed frequently when using PAML programs. As of the current version of
    PAML, the chi2 program does not allow passing both a test statistic and the
    degrees of freedom as command-line arguments.

    """
#    df = 2
#    statistic = 7.21
#    cdf_chi2(df, statistic)

