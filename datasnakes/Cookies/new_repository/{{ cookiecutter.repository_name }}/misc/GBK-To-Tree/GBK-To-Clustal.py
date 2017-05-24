# -*- coding: utf-8 -*-
"""

GBK To Clustal Version 1.0

Part A: Parse genbank files, extract the desired features, and store those
        features in fasta files or genbank files for downstream usage.

Part B: Align multifasta files using Clustal Omega and produce output
        in Phylip format.

Author: Shaurita D. Hutchins

Date created: January 31, 2017

*The current version of this program is designed for use on the MCSR
(a remote supercomputer).

"""

# List of modules used in this script
import sys
import time as t
import os
# This module is used for parsing genbank files & extracting record features.
from Bio import SeqIO
import csv  # Read a comma delimited list.
# The following is a command line wrapper for Clustal Omega
from Bio.Align.Applications import ClustalOmegaCommandline as ClustalCMD
#from Bio.Seq import Seq
#from Bio.Alphabet import generic_dna


# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------

# This prints a short description of the script.
print("#### GBK To Tree Version 2.0 ###" + (2 * "\n"))

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------


# Just a little fun. Hehe.
print("(•_•)" + "\n")
t.sleep(.75)
os.system("cls")
print("( •_•)>⌐■-■" + "\n")
t.sleep(.75)
os.system("cls")
print("(⌐■_■)")
t.sleep(.75)
print("\n" + "Let's GEAUX!!!!" + (2 * "\n"))
t.sleep(.75)
print("#" + (100 * "-") + "#" + "\n" + "#" + (100 * "-") + "#" + "\n" + "#" +
      (100 * "-") + "#" + "\n")
t.sleep(.75)

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------

# This is just a step to confirm that you are ready to start with Part A
# of this script.
print("In Part A, you will parse genbank files, extract the desired features," +
      " and store those features in fasta files or genbank files for downstream" +
      " usage.")
print("\n")

# Set Variable for home directory
home = '/ptmp/r2295/bin/Orthologs-Project/'
os.chdir(home)
# Set Variable for scripts dir
scripts = '/work5/r2295/bin/Orthologs-Project/Scripts/'

# Set Variable for scripts dir
docs = '/work5/r2295/bin/Orthologs-Project/ScriptFiles'
os.chdir(docs)

# Create and read a list of organisms from a .csv file.
name_list = []  # Initialize list of organisms
name_list.append('')
name = open('commonnames.csv')  # Open a comma delimited list of organisms.
file1 = csv.reader(name)
for cn in file1:  # Format a list of organisms
    cn = str(cn)
    cn = cn.replace("'", "")
    cn = cn.replace("[", "")
    cn = cn.replace("]", "")
    cn = cn.replace(" ", "_")
    name_list.append(cn)
print("List of Short Names:" + "\n")
print(name_list)  # Print the list of organisms

# -----------------------------------------------------------------------------

org_list = []  ## Initialize list of Organisms
org_list.append('')
o = open('Organisms.csv')
file2 = csv.reader(o)
for org in file2:    ##Format a list of organisms
    org = str(org)
    org = org.replace("'", "")
    org = org.replace("[", "")
    org = org.replace("]", "")
    org = org.replace(" ", "_")
    org_list.append(org)
print("List of Species Names:" + "\n")
print(org_list)  # Print the list of organisms
# -----------------------------------------------------------------------------

# Read a list of gene names by tier
tiers = ["Tier1", "Tier2", "Tier3", "Nontiered"]
files = ['f2', 'f3', 'f4', 'f5']
tier_count = 0
file_count = 0

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------

# First for loop to easily open and read through tier files.
# Creates tier directories.
for tier, file in zip(tiers, files):

    tier_count = tier_count + 1
    file_count = file_count + 1

    tierfile = open(tiers[0] + "genes.csv")
    file = csv.reader(tierfile)

    # Set a variable for the genbank file directory
    gbk_dir = '/ptmp/r2295/bin/Orthologs-Project/GBK-DIR/'
    a = gbk_dir  # Set variable for the genbank directories and files

    # Create the main CDS directory
    b = home + "CDS-DIR"
    # Create a directory or don't if it exists.
    os.makedirs('%s' % b, exist_ok=True)

    # Create Tier directories
    tdir = b + "/" + str(tiers[0]) + "-DIR"
    os.makedirs('%s' % tdir, exist_ok=True)


    # Let me know where I am right before I start the loop.
    print("\n" + "The current working directory is " +
          os.getcwd() + (2 * "\n"))  # Print current working directory
    Gene_count = 0

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------

    # 2nd loop. It creates gene directories for output.
    for Gene in file:
        Gene_count = Gene_count + 1

        # Create directories for different output files and designate variables as
        # paths.

        # Directories for gene-specific genbank files.
        d = a + "/" + str(tiers[0]) + "/" + str(Gene[0])

        # Create directories for CDS/Fasta files per gene
        c = tdir + "/" + str(Gene[0]) + "_CDS"  # Setting the variable for the path
        os.makedirs('%s' % c, exist_ok=True)

        # Create a directory for alignment files
        e = home + "/Alignments"
        os.makedirs('%s' % e, exist_ok=True)

        # Create directories for alignment files per gene
        # Setting the variable for the path
        f = e + "/" + str(Gene[0]) + "_Aligned"
        os.makedirs('%s' % f, exist_ok=True)

        # Set a variable for the PhyloAnalysis directory
        phylo_dir = '/ptmp/r2295/bin/Orthologs-Project/PhyloAnalysis/'
        g = phylo_dir  # Tag a shorter variable as the PhyloAnalysis directory path

        # Create directories for PhyMl output
        h = g + str(Gene[0]) + "_PhyloTrees"
        os.makedirs('%s' % h, exist_ok=True)

        # Create directories for PAML files
        i = g + str(Gene[0]) + "_PAML_Output"
        os.makedirs('%s' % i, exist_ok=True)

        # Change to a genbank file directory of the gene (will loop through list
        # of genes specified)
        os.chdir(d)
        os.listdir(d)  # Make a list of the files in the current directory

        # Print current working directory
        print("➜ Current gene directory: " + os.getcwd() + "\n")
        # input("    If this is the desired directory, press ENTER.")
        print("\n")

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------

        # Part A: Parse genbank files, extract the desired features, and store
        # those features in fasta files or genbank files for downstream usage.

        # Loop that establishes organism list, reads genbank file, and
        # creates/opens new fasta file.
        file_count = 0
        for Organism, Name in zip(org_list, name_list):
            file_count = file_count + 1
            maximum = 0
            if Organism == '':
                continue
            os.chdir(d)  # Directory of genbank files
            record = SeqIO.read(str(Gene[0]) + "_" + Organism + ".gbk", "genbank")
            os.chdir(c)  # Change to directory for cds.fasta and aa.fasta files
            # You can also create a genbank file.
            output_handle = open(str(Gene[0]) + "_" + Organism + "_cds_nucl.fasta", "w")
            count = 0

# -----------------------------------------------------------------------------

            # Loop that extracts specific features and writes them to previously
            # created file.
            for feature in record.features:
                # Other annotated features are 'Gene', 'mRNA', 'CDS', and 'ncRNA'.
                if feature.type == "CDS":
                    count = count + 1
                    # Use record.dbxrefs here. Look up record features in Ipython
                    # using 'dir(record)'.
                    feature_name = Name
                    feature_seq = feature.extract(record.seq)
                    # Simple FASTA output without line wrapping:
                    output_handle.write(
                        ">" + feature_name + "  " + "\n" + str(feature_seq) + "\n")
                    output_handle.close()
                    print(Organism + "\n" + feature_name + "\n" + str(count) +
                          " CDS sequence was extracted from " + Organism + "." + (2 * "\n"))


#                    # Translate the sequence to an amino acid sequence as well
#                    coding_dna = Seq(str(record.seq), generic_dna)
#                    # Translate the nucleotide sequence & save to a separate file
#                    with open(str(Gene[0]) + "_" + Organism + "_cds_aa.fasta") as aa:
#                        aa.write(str(">" + feature_name + "\n" +
#                                     coding_dna.translate(table=1, to_stop=True)))

# -----------------------------------------------------------------------------
            # Creating space between output
            print("#" + (100 * "-") + "#" + "\n" + "#" + (100 * "-") + "#" + "\n")
            print("\n")
        print((50 * "✓") + "\n")
        print("\n")

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------

        # Part B: Align multifasta files using Clustal Omega and produce output in
        # Phylip format.

        # Move to CDS directory
        os.chdir(c)  # Change to the CDS files directory
        os.listdir(c)  # Make a list of the files in the current directory

        # Print current working directory
        print("➜ Current CDS/Gene directory: " + os.getcwd() + "\n")
        # input("If this is the desired directory, press ENTER.")
        print("\n")

        # Echos all commands in the current shell.
        os.system("set -x")

        # Creates a copy of the Homo sapiens cds file and renames it as a profile
        # sequence.
        os.system("cp " + str(Gene[0]) + "_Homo_sapiens_cds_nucl.fasta profile_nucl.fasta")
        print("\n")

        # Remove the Homo sapiens cds.fasta file before concatenation.
        os.remove(str(Gene[0]) + "_Homo_sapiens_cds_nucl.fasta")
        print("\n")

#        # Remove the Homo sapiens cds.fasta file before concatenation.
#        os.remove(str(Gene[0]) + "_Homo_sapiens_cds_aa.fasta")
#        print("\n")

        # Uses command line to concatenate nucleotide fasta files in current directory.
        os.system("cat *_cds_nucl.fasta* > " + str(Gene[0]) + "_cds_nucl.fasta")
        print("\n")

#        # Uses command line to concatenate amino acid fasta files in current directory.
#        os.system("cat *_cds_aa.fasta* > " + str(Gene[0]) + "_cds_aa.fasta")
#        print("\n")

        # Copies the profile.fasta and concatenated cds.fasta file to output dir.
        os.system("cp {" + str(Gene[0]) + "_cds_nucl.fasta,profile_nucl.fasta} " + c + "/")

#        # Copies the profile.fasta and concatenated cds.fasta file to output dir.
#        os.system("cp {" + str(Gene[0]) + "_cds_aa.fasta,profile_aa.fasta} " + c + "/")

        # Directory change to output directory
        os.chdir(f)

        # Run Clustal Omega and produce output in phylip format
        print("\n")
        print("Clustal Ω will align the sequences and produce output in phylip format.")
        print("\n")

# -----------------------------------------------------------------------------

        # Align the nuclueotide sequences
        in_file1 = c + "/" + str(Gene[0]) + "_cds_nucl.fasta"
        out_file1 = str(Gene[0]) + "_cds_nucl_aligned"
        clustalo_cline1 = ClustalCMD(profile1= c + "/profile_nucl.fasta", infile=in_file1,
                                     outfile=out_file1, seqtype="DNA",
                                     infmt="fasta", outfmt="phy", iterations=8,
                                     distmat_full_iter=True,
                                     verbose=True,
                                     threads=8, force=True,
                                     log=str(Gene[0]) + "_nucl_fasta_log.txt")
        stdout1, stderr1 = clustalo_cline1()

# -----------------------------------------------------------------------------

#        # Align the amino acid sequences
#        in_file2 = str(Gene[0]) + "_cds_aa.fasta"
#        out_file2 = str(Gene[0]) + "_cds_aa_aligned"
#        clustalo_cline2 = ClustalCMD(profile1="profile_aa.fasta", infile=in_file2,
#                                     outfile=out_file2, seqtype="protein",
#                                     infmt="fasta", outfmt="fasta, phy", iterations=8,
#                                     distmat_full_iter=True,
#                                     verbose=True,
#                                     threads=8, force=True,
#                                     log=str(Gene[0]) + "_aa_fasta_log.txt")
#        stdout2, stderr2 = clustalo_cline2()

# -----------------------------------------------------------------------------
        # Prints the standard output and error of the clustal omega command
        clustalo_cline1()
        print(stdout1, stderr1)
        print("Fasta formatted nucleotide alignment file has been created for  " + str(Gene[0]) + "." + "\n")
#        clustalo_cline2()
#        print(stdout2, stderr2)
#        print("Fasta formatted amino acid alignment file has been created for " + str(Gene[0]) + "." + "\n")


        # Change to cds.fasta file directory of current gene
        os.chdir(c)

        # Create a copy of the profile_nucl.fasta file and rename it to original name.
        os.system("cp profile_nucl.fasta " + str(Gene[0]) + "_Homo_sapiens_cds_nucl.fasta")

#        # Create a copy of the profile_aa.fasta file and rename it to original name.
#        os.system("cp profil_aa.fasta " + str(Gene[0]) + "_Homo_sapiens_cds_aa.fasta")


        # Remove the profile.fasta file after creating alignments
        os.remove("profile_nucl.fasta")
#        os.remove("profile_aa.fasta")



print("This script has finished. ✓ You may exit now.")
o.close()
g.close()
sys.exit()
