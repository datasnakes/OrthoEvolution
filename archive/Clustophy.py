# -*- coding: utf-8 -*-
"""

@author: S. Hutchins
"""
# Align multiple sequences located in fasta files using alignment
# via MCSR and create a .phylip output format.

import csv
# List of modules used
import os
import sys
import time

from Bio.Align.Applications import ClustalOmegaCommandline

# Mark start of script with printed text.
print("\n" + (70 * "#") + "\n" + "#### Align fasta files using alignment via the MCSR.  ####" + "\n" + (70 * "#") + "\n")

# Read a list of gene names from a .csv file.
g = open('Gene_1_names.csv')  # 1st column - gene names
file1 = csv.reader(g)

# Let me know where I am right before I start the loop.
print("\n" + "The current working directory is "+ os.getcwd() + (2 * "\n"))  # Print current working directory
Gene_count = 0

from Orthologs.manager.comparative_genetic_analysis import CompGenAnalysis
x = CompGenAnalysis.get_master_lists()
for Gene in file1:
    Gene_count = Gene_count + 1

    # Assess current working directories and establishing home and output directories.

    a = '/ptmp/r2295/bin/CDS/' # Home directory
    home = a  # Location of CDS.fasta directories and files
    b = '/ptmp/r2295/bin/CDS/Alignments/' # Output directory
    output = b  #Location of aligned alignments files
    os.chdir(output) # Directory change to output directory

    # Create directories for alignment files.
    c = b + "./" + str(Gene[0]) + "_Aligned"
    os.makedirs('%s' % c, exist_ok=True) # Create a directory or don't if it exists.

    # Create directories for CDS files if they don't exist.
    os.chdir(home)  # Change to the home directory
    d = a + "./" + str(Gene[0])
    os.makedirs('%s' % d, exist_ok=True)

    os.chdir(d) # Change to cds.fasta file directory
    os.listdir(d) # Make a list of the files in the current directory
    print("➜ Current CDS/Gene directory: "+ os.getcwd() + "\n")  # Print current working directory
    time.sleep(.3)
    #input("    If this is the desired directory, press ENTER.")
    print("\n")

    # Echos all commands in the current shell.
    os.system("set -x")

    # Creates a copy of the Homo sapiens cds file and renames it as a profile sequence.
    os.system("cp " + str(Gene[0]) + "_Homo_sapiens_cds.fasta profile.fasta")
    print("\n")

    # Remove the Homo sapiens cds.fasta file before concatenation.
    os.remove(str(Gene[0]) + "_Homo_sapiens_cds.fasta")
    print("\n")

    # Uses command line to concatenate fasta files in current directory.
    os.system("cat *_cds.fasta* > " + str(Gene[0]) + "_cds.fasta")
    print("\n")

    # Uses command line to remove "PREDICTED: " from beginning of lines in fasta file.
    os.system("sed -i 's/PREDICTED: //g' " +  str(Gene[0]) + "_cds.fasta")
    print("\n")

    # Copies the profile.fasta and concatenated cds.fasta file to output dir.
    os.system("cp {" + str(Gene[0]) + "_cds.fasta,profile.fasta} " + c + "/")

    # Directory change to output directory
    os.chdir(c)

    # Output in phylip format
    print("\n" + "Clustal Ω will align the sequences and produce output in phylip format." + "\n")
    in_file1 = str(Gene[0]) + "_cds.fasta"
    out_file1 = str(Gene[0]) + "_aligned.phy"
    clustalo_cline1 = ClustalOmegaCommandline(profile1="profile.fasta", infile=in_file1, outfile=out_file1, seqtype="DNA",
                                             infmt="fasta", outfmt="phylip", iterations=4, distmat_full_iter=True, verbose=True,
                                             threads=8, force=True, log=str(Gene[0]) + "_phy_log.txt")
    stdout1, stderr1 = clustalo_cline1()
    clustalo_cline1()
    print(stdout1, stderr1)
    print("Phylip formatted alignment file has been created." + "\n")

    # Change to cds.fasta file directory of current gene
    os.chdir(d)

    # Create a copy of the profile.fasta file and rename it to original name.
    os.system("cp profile.fasta " + str(Gene[0]) + "_Homo_sapiens_cds.fasta")

    # Remove the profile.fasta file after creating alignments
    os.remove("profile.fasta")

sys.exit("✓✓✓✓ This script has completed. ✓✓✓✓")  # Exit the script.
