# -*- coding: utf-8 -*-
"""
Created on Fri Sep 16 11:31:04 2016

@author: Shaurita D. Hutchins
"""
import time
print("This program is designed to parse genbank files and extract desired features.")
time.sleep(1.5)
input("If you would like to start this program, press enter." )

# Listing of modules that will be used in loop(s).
#import time
import os
from Bio import SeqIO #Downstream usage of this module for parsing genbank files & extracting record features.
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
import csv # Read list of files.

# Assessing current working directories and establishing home and output directories.
print("The current working directory is "+ os.getcwd() + "\n")  # Print current working directory
time.sleep(0.5)

x = r'C:\Users\shutchins2\Desktop\In Progress\Code\GBK2TREE'
os.chdir(x)

org_list = []  ## Initialize list of Organisms
org_list.append('')
o = open('Organisms.csv')
file1 = csv.reader(o)
for org in file1:    ##Format a list of organisms
    org = str(org)
    org = org.replace("'", "")
    org = org.replace("[", "")
    org = org.replace("]", "")
    org = org.replace(" ", "_")
    org_list.append(org)
print(org_list)
time.sleep(0.5)

g = open('Gene_1_names.csv')  ## 1st column - gene names
file2 = csv.reader(g)

Gene_count = 0

for Gene in file2:
    Gene_count = Gene_count + 1


    a = r'C:/Users/shutchins2/Desktop/In Progress/Code/GBK2TREE/GBK' #Home Directory
    home = a    #Location of genbank directories and files
    b = r'C:/Users/shutchins2/Desktop/In Progress/Code/GBK2TREE/CDS' # Output directory of cds.fasta files
    output = b
    os.chdir(output) ## Directory Change: Output directory

    # Create/access directories for CDS/Fasta files.
    c = b + "./" + str(Gene[0])
    os.makedirs('%s' % c, exist_ok=True) # Create directory or don't if it exists.

    # Create/access directories for Genbank files.
    os.chdir(home) ## Make a list of files
    d = a + "./" + str(Gene[0])
    os.makedirs('%s' % d, exist_ok=True)

    os.chdir(d) # Change to genbank file directory
    os.listdir() # Make a list of the files in the current directory


##############################################################################################
##############################################################################################
##############################################################################################
##############################################################################################


##### Part 1: Parse genbank files & write CDS to 1 fasta file per gene for multiple alignments. #####
    print("##### Parse genbank files & write CDS or other features to 1 fasta file per gene for multiple alignments. #####")
    time.sleep(.3)
    file_count = 0
    for Organism in org_list: # GBK to CDS (or other feature) loop
        file_count = file_count + 1
        maximum = 0
        if Organism == '':
            continue
        os.chdir(d) # Directory of genbank files
        record = SeqIO.read(str(Gene[0]) + "_" + Organism + ".gbk", "genbank")
        os.chdir(c) #Change to directory for cds.fasta files
        output_handle = open(str(Gene[0]) + "_" + Organism + "_cds.fasta", "w")
        count = 0
        for feature in record.features:
        # Other annotated features are 'Gene', 'mRNA', 'CDS', and 'ncRNA'.
            if feature.type == "CDS":
                count = count + 1
                # Use record.dbxrefs here. Look up record features in Ipython
                # using 'dir(record)'.
                feature_name = Organism
                feature_seq = feature.extract(record.seq)

                # Simple FASTA output without line wrapping:
                output_handle.write(
                    ">" + feature_name + "  " + "\n" + str(feature_seq) + "\n")
                output_handle.close()
                print(Organism + "\n" + feature_name + "\n" + feature_seq + "\n" + "\n" + str(
                    count) + " CDS sequence was extracted from " + Organism + "." + (2 * "\n"))

#                # Translate the sequence to an amino acid sequence as well
#                coding_dna = Seq(record.seq, generic_dna)
#                # Translate the nucleotide sequence & save to a separate file
#                with open(str(Gene[0]) + "_" + Organism + "_cds_aa.fasta") as aa:
#                    aa.write(str(">" + feature_name + "\n" +
#                                 coding_dna.translate(table=1, to_stop=True)))

print("This script has finished. You may exit now.")
o.close()
g.close()
exit()