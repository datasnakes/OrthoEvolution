# -*- coding: utf-8 -*-
"""

GBK To Phylip Version 1.0

Part A: Parse genbank files, extract the desired features, and store those
        features in fasta files or genbank files for downstream usage.

Part B: Align multifasta files using Clustal Omega and produce output
        in Phylip format.

Part C: Use Phylip to produce phylogenetic trees.


Author: Shaurita D. Hutchins

Date created: January 19, 2017


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
import pexpect
#import logging as log
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio.Phylo.Applications import PhymlCommandline

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------

# This prints a short description of the script.
print("#### GBK To Tree Version 2.0 ###" + "\n")

# There are a few short "sleeps" in this program so that it can be completed
# interactively and read in real time.
t.sleep(.5)

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
print((107 * "#") + "\n" + (107 * "#") + "\n" + (107 * "#") + (3 * "\n"))
t.sleep(.75)

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------

# This is just a step to confirm that you are ready to start with Part A
# of this script.
input("In Part A, you will parse genbank files, extract the desired features," +
      " and store those features in fasta files or genbank files for downstream" +
      " usage. If you would like to start Part A, press enter. ")
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
t.sleep(0.5)

# -----------------------------------------------------------------------------

# Read a list of gene names by tier
tiers = ["tier1", "tier2", "tier3", "nontiered"]
files = [2, 3, 4]

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------

# First for loop to easily open and read through tier files.
for tier, file in zip(tiers, files):
    tier[0] = open(tier[0] + "genes.csv")
    file[0] = csv.read(tier[0])


    # Let me know where I am right before I start the loop.
    print("\n" + "The current working directory is " +
          os.getcwd() + (2 * "\n"))  # Print current working directory
    t.sleep(.5)
    Gene_count = 0

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------

    # 2nd loop. It creates gene directories for output.
    for Gene in file[0]:
        Gene_count = Gene_count + 1

        # Create directories for different output files and designate variables as
        # paths.

        # Set Variable for home directory
        home = '/ptmp/r2295/bin/Orthologs-Project/'

        # Set a variable for the genbank file directory
        gbk_dir = '/ptmp/r2295/bin/Orthologs-Project/GBK-DIR/'
        a = gbk_dir  # Set variable for the genbank directories and files

        # Create the main CDS directory
        b = home + "./CDS-DIR"
        # Create a directory or don't if it exists.
        os.makedirs('%s' % b, exist_ok=True)

        # Create directories for CDS/Fasta files per gene
        c = b + "./" + str(Gene[0]) + "_CDS"  # Setting the variable for the path
        os.makedirs('%s' % c, exist_ok=True)

    #    # Create directories for Amino Acid fasta files per gene
    #    aa = b + "./" + str(Gene[0]) + "_AminoAcid"
    #    os.makedirs('%s' % c, exist_ok=True)

        # Create a directory for alignment files
        e = home + "./Alignments"
        os.makedirs('%s' & e, exist_ok=True)

        # Create directories for alignment files per gene
        # Setting the variable for the path
        f = e + "./" + str(Gene[0]) + "_Aligned"
        os.makedirs('%s' % f, exist_ok=True)

        # Set a variable for the PhyloAnalysis directory
        phylo_dir = '/ptmp/r2295/bin/Orthologs-Project/PhyloAnalysis/'
        g = phylo_dir  # Tag a shorter variable as the PhyloAnalysis directory path

        # Create directories for PhyMl output
        h = g + "./" + str(Gene[0]) + "_PhyloTrees"
        os.makedirs('%s' % h, exist_ok=True)

        # Create directories for PAML files
        i = g + "./" + str(Gene[0]) + "_PAML_Output"
        os.makedirs('%s' % i, exist_ok=True)

        # Change to the genbank file directory to begin Part A
        # Set a variable to the genbank file directory per gene - subdirectory
        d = a + "./" + str(Gene[0])
        # Change to a genbank file directory of the gene (will loop through list
        # of genes specified)
        os.chdir(d)
        os.listdir(d)  # Make a list of the files in the current directory

        # Print current working directory
        print("➜ Current gene directory: " + os.getcwd() + "\n")
        t.sleep(.3)
        # input("    If this is the desired directory, press ENTER.")
        print("\n")

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------

        # Part A: Parse genbank files, extract the desired features, and store
        # those features in fasta files or genbank files for downstream usage.

        t.sleep(.3)
        file_count = 0

# -----------------------------------------------------------------------------

        # Loop that establishes organism list, reads genbank file, and
        # creates/opens new fasta file.
        for Organism in org_list:
            file_count = file_count + 1
            maximum = 0
            if Organism == '':
                continue
            os.chdir(d)  # Directory of genbank files
            record = SeqIO.read(str(Gene[0]) + "_" + Organism + ".gbk", "genbank")
            os.chdir(c)  # Change to directory for cds.fasta files
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
                    feature_name = (Organism)
                    feature_seq = feature.extract(record.seq)
                    # Simple FASTA output without line wrapping:
                    output_handle.write(
                        ">" + feature_name + "\n" + str(feature_seq) + "\n")
                    output_handle.close()
                    print(Organism + "\n" + feature_name + "\n" + feature_seq + "\n" + "\n" + str(
                        count) + " CDS sequence was extracted from " + Organism + "." + (2 * "\n"))
                    # t.sleep(0.15)

                    # Translate the sequence to an amino acid sequence as well
                    coding_dna = Seq(record.seq, generic_dna)
                    # Translate the nucleotide sequence & save to a separate file
                    with open(str(Gene[0]) + "_" + Organism + "_cds_aa.fasta") as aa:
                        aa.write(str(">" + feature_name + "\n" +
                                     coding_dna.translate(table=1, to_stop=True)))

# -----------------------------------------------------------------------------

            print((100 * "#") + "\n" + (100 * "#") + "\n" +
                  (100 * "#") + "\n")  # Creating space between output
            print(2 * "\n")
        print((50 * "★") + "\n")
        # This input lets you know which gene is next.
        # input("If you would like to continue to the next gene, press ENTER. ")
        print("\n")

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------

        # Part B: Align multifasta files using Clustal Omega and produce output in
        # Phylip format.

        # Move to CDS directory
        os.chdir(c)  # Change to the CDS files directory
        os.listdir(c)  # Make a list of the files in the current directory
        t.sleep(.3)

        # Print current working directory
        print("➜ Current CDS/Gene directory: " + os.getcwd() + "\n")
        # input("If this is the desired directory, press ENTER.")
        print("\n")

        # Echos all commands in the current shell.
        os.system("set -x")

        # Creates a copy of the Homo sapiens cds file and renames it as a profile
        # sequence.
        os.system("cp " + str(Gene[0]) + "_Homo_sapiens_cds.fasta profile_nucl.fasta")
        print("\n")

        # Remove the Homo sapiens cds.fasta file before concatenation.
        os.remove(str(Gene[0]) + "_Homo_sapiens_cds_nucl.fasta")
        print("\n")

        # Remove the Homo sapiens cds.fasta file before concatenation.
        os.remove(str(Gene[0]) + "_Homo_sapiens_cds_nucl.fasta")
        print("\n")

        # Uses command line to concatenate nucleotide fasta files in current directory.
        os.system("cat *_cds_nucl.fasta* > " + str(Gene[0]) + "_cds_nucl.fasta")
        print("\n")

        # Uses command line to concatenate amino acid fasta files in current directory.
        os.system("cat *_cds_aa.fasta* > " + str(Gene[0]) + "_cds_aa.fasta")
        print("\n")

        # Copies the profile.fasta and concatenated cds.fasta file to output dir.
        os.system("cp {" + str(Gene[0]) + "_cds_nucl.fasta,profile_nucl.fasta} " + c + "/")

        # Copies the profile.fasta and concatenated cds.fasta file to output dir.
        os.system("cp {" + str(Gene[0]) + "_cds_aa.fasta,profile_aa.fasta} " + c + "/")

        # Directory change to output directory
        os.chdir(f)

        # Run Clustal Omega and produce output in phylip format
        print("\n")
        print("Clustal Ω will align the sequences and produce output in phylip format.")
        print("\n")

# -----------------------------------------------------------------------------

        # Align the nuclueotide sequences
        in_file1 = str(Gene[0]) + "_cds_nucl.fasta"
        out_file1 = str(Gene[0]) + "_cds_nucl_aligned.fasta"
        clustalo_cline1 = ClustalCMD(profile1="profile_nucl.fasta", infile=in_file1,
                                     outfile=out_file1, seqtype="DNA",
                                     infmt="fasta", outfmt="fasta", iterations=4,
                                     distmat_full_iter=True,
                                     verbose=True,
                                     threads=8, force=True,
                                     log=str(Gene[0]) + "_nucl_fasta_log.txt")
        stdout1, stderr1 = clustalo_cline1()

# -----------------------------------------------------------------------------

        # Align the amino acid sequences
        in_file2 = str(Gene[0]) + "_cds_aa.fasta"
        out_file2 = str(Gene[0]) + "_cds_aa_aligned.fasta"
        clustalo_cline2 = ClustalCMD(profile1="profile_aa.fasta", infile=in_file2,
                                     outfile=out_file2, seqtype="protein",
                                     infmt="fasta", outfmt="fasta", iterations=4,
                                     distmat_full_iter=True,
                                     verbose=True,
                                     threads=8, force=True,
                                     log=str(Gene[0]) + "_aa_fasta_log.txt")
        stdout2, stderr2 = clustalo_cline2()

# -----------------------------------------------------------------------------
        # Prints the standard output and error of the clustal omega command
        clustalo_cline1()
        print(stdout1, stderr1)
        print("Fasta formatted nucleotide alignment file has been created for  " + str(Gene[0]) + "." + "\n")
        clustalo_cline2()
        print(stdout2, stderr2)
        print("Fasta formatted amino acid alignment file has been created for " + str(Gene[0]) + "." + "\n")


        # Change to cds.fasta file directory of current gene
        os.chdir(c)

        # Create a copy of the profile_nucl.fasta file and rename it to original name.
        os.system("cp profile_nucl.fasta " + str(Gene[0]) + "_Homo_sapiens_cds_nucl.fasta")
        # Create a copy of the profile_aa.fasta file and rename it to original name.
        os.system("cp profil_aa.fasta " + str(Gene[0]) + "_Homo_sapiens_cds_aa.fasta")


        # Remove the profile.fasta file after creating alignments
        os.remove("profile_nucl.fasta"), os.remove("profile_aa.fasta")

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------

        # Part C: Use Phylip to produce phylogenetic trees

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

