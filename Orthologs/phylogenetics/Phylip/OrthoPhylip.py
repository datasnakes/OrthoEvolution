# -*- coding: utf-8 -*-
"""
Created on Fri Oct 28 2016

@author: S. Hutchins

Note: I'm certain that there's a way to truncate this script.
"""
# Uses multiple sequence alignments in phylip format to create
# phylogenetic trees and other phylogenetic analysis.

# List of modules
import csv
import os
import sys
import pexpect  # I used this to feed input into shell executable
import logging as log

# Mark start of program with printed text description/title.
print("\n" + (81 * "#") + "\n")
print("#### Use multiple sequence alignments to create phylogenetic trees using Phylip. ####")
print("\n" + (81 * "#") + "\n")


# Create a log file for this script
FORMAT = '%(asctime)s | %(levelname)s: %(message)s'
log.basicConfig(
    filename='OrthoPhylip.log',
    format=FORMAT,
    datefmt='Date: %m/%d/%Y Time: %I:%M:%S %p',
    level=log.INFO)

# Read a list of gene names from a .csv file.
genes = open('Gene_1_names.csv')  # 1st column - List of genes used
file1 = csv.reader(genes)
log.info('The script is beginning.')
log.info('This is the list of genes: ' + print(file1))

# Let me know where I am right before I start the loop.
print("\n" + "The current working directory is " +
      os.getcwd() + (2 * "\n"))  # Print current working directory
Gene_count = 0

# -----------------------------------------------------------------------------

for Gene in file1:  # Loop to create trees for each gene related aligned file
    Gene_count = Gene_count + 1

    # Assess current working directories and establishing home and output
    # directories.
    a = '/ptmp/r2295/bin/Genes_1_CDS/Alignments/'  # Home directory
    home = a  # Location of alignment directories and files
    b = '/ptmp/r2295/bin/Genes_1_CDS/Alignments/PhyloAnalysis/'  # Output dir
    output = b  # Location of phylogenetic analyses files
    os.chdir(output)  # Directory change to output directory

    # Create directories for phylogenetic tree & analyses files.
    c = b + "./" + str(Gene[0]) + "_PhyloTrees"
    # Create a directory or don't if it exists.
    os.makedirs('%s' % c, exist_ok=True)

    # Create directories for alignment files if they don't exist.
    os.chdir(home)  # Change to the home directory
    d = a + "./" + str(Gene[0]) + "_Aligned"
    os.makedirs('%s' % d, exist_ok=True)

    os.chdir(d)  # Change to alignment file directory
    os.listdir(d)  # Make a list of the files in the current directory
    # Print current working directory
    print("➜ Current CDS/Gene directory: " + os.getcwd() + "\n")
    input("If this is the desired directory, press ENTER.")
    print("\n")

# -----------------------------------------------------------------------------

    # Echos all commands in the current shell.
    os.system("set -x")

    # Creates a copy of the phylip alignnment file and renames it as an infile.
    os.system("cp " + str(Gene[0]) + "_aligned.phy infile")
    print("\n")

    # Copies infile to the output directory
    os.system("cp infile " + c + "/")

    # Directory change to output directory
    os.chdir(c)

    # Create a variable for os.rename
    rn = os.rename

# -----------------------------------------------------------------------------
    # Maximum Likelihood using Phylip executable, dnaml, within unix shell
    dnaml = pexpect.spawnu("dnaml infile")
    dnaml.sendline("Y\r")
    dnaml.waitnoecho()
    rn('"outfile, "' + str(Gene[0]) + '_maxlike"')
    rn('"outtree, "' + str(Gene[0]) + '_maxliketree"')

# -----------------------------------------------------------------------------
    # Maximum Parsimony using Phylip executable, dnapars, within unix shell
    dnapars = pexpect.spawnu("dnapars infile")
    dnapars.sendline("Y\r")
    dnapars.waitnoecho()
    rn('"outfile, "' + str(Gene[0]) + '_maxpars"')
    rn('"outtree, "' + str(Gene[0]) + '_maxparstree"')

# -----------------------------------------------------------------------------
    # Distance Matrix using the Phylip executable, dnadist, within unix shell
    dnadist = pexpect.spawnu("dnadist infile")
    dnadist.sendline("Y\r")
    dnadist.waitnoecho()
    rn('"outfile", "' + str(Gene[0]) + '_dnadist"')

# -----------------------------------------------------------------------------

log.info('This script has finished. Check your phylogenetic trees.')
sys.exit("This script has finished. Check your phylogenetic trees.")
