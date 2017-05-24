# -*- coding: utf-8 -*-
"""
Created on Mon Jan  9 13:15:22 2017

@author: S. Hutchins

Simple test for clustal omega to phyml command line to PAML on the MCSR
"""
# List of modules used
import os
import sys
import csv
from Bio.Align.Applications import ClustalOmegaCommandline
from Bio.Phylo.Applications import PhymlCommandline
from Bio import AlignIO
from Bio import Phylo

home = '/work5/r2295/bin/Orthologs-Project'
os.chdir(home)
# This is just a step to confirm that you are ready to start.
input("This script is designed to work on Windows. If you would like to start this script, press enter. ")
print("\n")

# Create and read a list of organisms from a .csv file.
org_list = []  # Initialize list of organisms
org_list.append('')
o = open('Organisms.csv')  # Open a comma delimited list of organisms.
file1 = csv.reader(o)
for org in file1:  # Format a list of organisms
    org = str(org)
    org = org.replace("'", "")
    org = org.replace("[", "")
    org = org.replace("]", "")
    org = org.replace(" ", "_")
    org_list.append(org)
print("List of Organisms" + "\n")
print(org_list)  # Print the list of organisms


g = open('tier1genes.csv')  # 1st column - gene names
file2 = csv.reader(g)

# Let me know where I am right before I start the loop.
print("\n" + "The current working directory is " +
      os.getcwd() + (2 * "\n"))  # Print current working directory

Gene_count = 0

# Loop to create gene directories. This loop includes 2 other loops. Be
# mindful.
for Gene in file2:
    Gene_count = Gene_count + 1

    # Directories
    clustalo = r'C:\Users\shutchins2\Desktop\Executables\clustal-omega-1.2.2-win64'
    phyml_dir = r''
    cds_dir = r'C:\Users\shutchins2\Desktop\In Progress\Code\GBK2TREE\CDS' + "\\" +  Gene[0] + "\\"
    cdsmain = r'C:\Users\shutchins2\Desktop\In Progress\Code\GBK2TREE\CDS'

    # Change to home directory
    os.chdir(home)

    # Creates a temporary directory for the profile sequence/hmm and move it to that directory.
    os.system("mkdir test_output")
    test_out = home + "\\test_output\\" + Gene[0]
    os.makedirs('%s' % test_out, exist_ok=True)
    os.chdir(cds_dir)
    print(os.getcwd())
    #(print(os.listdir(cds_dir)))
    #input("Continue?" )
    x = r'C:\Users\shutchins2\Desktop\"In Progress"\Code\GBK2TREE\"Alignments & Trees"\"Alignment Scripts and Files"\test_output' + "\\" +  Gene[0] + "\\"
    os.system("move " + Gene[0] + "_Homo_sapiens_cds.fasta " + x)


    # Uses command line to concatenate fasta files in current directory.
    os.system("copy *_cds.fasta* " + Gene[0] + "_cds.fasta")
    z = r'C:\Users\shutchins2\Desktop\"In Progress"\Code\GBK2TREE\CDS'
    os.system("move " + Gene[0] + "_cds.fasta " + z)
    os.chdir(cdsmain)

    #------------------------------------------------------------------------------
    # Clustal Omega
    #------------------------------------------------------------------------------

    # Mark start of program with printed text.
    print("\n" + (70 * "#") + "\n" + "#### Align fasta files using Clustal Omega.  ####" + "\n" + (70 * "#") + "\n")

    #input("If you would like to create a fasta file, press ENTER. ")
    # Output in fasta format
    print("\n" + "Clustal Ω will align the sequences and produce output in fasta format." + "\n")
    in_file3 = Gene[0] + "_cds.fasta"
    out_file3 = Gene[0] + "_cds_aligned.fasta"
    clustalo_cline3 = ClustalOmegaCommandline(profile1=test_out + "\\" + Gene[0] + "_Homo_sapiens_cds.fasta",
                                              infile=in_file3, outfile=out_file3, seqtype="DNA",
                                             infmt="fasta", outfmt="fasta", iterations=1,
                                             distmat_full_iter=True, verbose=True,
                                             force=True, log=Gene[0] + "_fasta_log.txt")

    stdout3, stderr3 = clustalo_cline3()
    clustalo_cline3()
    print(stdout3, stderr3)
    print("Fasta formatted alignment file has been created." + "\n")

    os.system("move " + Gene[0] + "_cds_aligned.fasta " + x)

    #------------------------------------------------------------------------------
    # PhyML
    #------------------------------------------------------------------------------
    # Mark start of program with printed text.
    print("\n" + (70 * "#") + "\n" + "#### Create phylogenetic trees using PhyML.  ####" + "\n" + (70 * "#") + "\n")


    os.chdir(home)

    # Use the phyml executable file
    phyml_exe = None
    exe_name = "PhyML-3.1_win32.exe" if sys.platform == "win32" else "phyml"
    phyml_exe = exe_name

    os.chdir(test_out)

    # Convert the file to relaxed-phylip format
    print("Convert the file to relaxed-phylip format.")
    AlignIO.convert(Gene[0] + "_cds_aligned.fasta", "fasta", Gene[0] + "_aligned.phy", "phylip-relaxed")
    y = r'C:\Users\shutchins2\Desktop\"In Progress"\Code\GBK2TREE\"Alignments & Trees"\"Alignment Scripts and Files"' + "\\"
    os.system("move " + Gene[0] + "_aligned.phy " + y)
    print("The file has been converted to relaxed-phylip format.")

    os.chdir(home)
    # Create the command & run phyml
    # Input a phylip formatted alignment file and describe the datatype ('nt' or 'aa')
    run_phyml = PhymlCommandline(phyml_exe, input=Gene[0] + '_aligned.phy', datatype='nt')
    print("\n" + "The following PhyML command is being run: "), print(run_phyml)
    out_log, err_log = run_phyml()


    #tree = Phylo.read('HTR1D_aligned.phy', 'newick')
    #Phylo.draw_ascii(tree)

    # Move phyml output files to appropriate directory
    os.system("move " + Gene[0] + "_aligned._phyml_tree.txt " + x)
    os.system("move " + Gene[0] + "_aligned._phyml_stats.txt " + x)
    os.system("move " + Gene[0] + "_aligned.phy " + x)
    print("PhyML tree files and stats have been created." + "\n")


    #------------------------------------------------------------------------------
    # PAML
    #------------------------------------------------------------------------------

#    # Change to the PAML executables directory
#    PAML = r'C:\Users\shutchins2\Desktop\Executables\paml4.9c'
#    os.chdir(PAML)
#
#    # Use the command line to call/run the executable programs
#    codeml = 'bin/codeml'
#    baseml = 'bin/baseml'
#    pamp = 'bin/pamp'
#
#    # Create a file for the standard output to be stored in
#    with open('codeml-results.txt', 'a') as codeml-results:
#        codeml-results.write(str(cmd(codeml)))
#
#    # Move the results file to the output directoy
#    os.system('cp codeml-results.txt ' + paml_dir + 'codeml-results.txt')
#    os.remove('codeml-log.txt')




