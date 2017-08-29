# -*- coding: utf-8 -*-
"""This script is designed to create a gi list based on the refseq_rna database
for each taxonomy id on the MCSR. It will also convert the gi list into a
binary file which is more efficient to use with NCBI's Standalone Blast tools.
"""

#------------------------------------------------------------------------------
# List of modules used.
#------------------------------------------------------------------------------

import os
import pandas as pd
import multiprocessing


#------------------------------------------------------------------------------
# Mark start of program with printed text.
#------------------------------------------------------------------------------

print("\n" + (90 * "#") + "\n" + "Create a GI list for each organism using the \
taxonomy id and the blastdbcmd tool on the MCSR.  ####" + "\n" + (90 * "#") + "\n")

#------------------------------------------------------------------------------
# Establish home directory
#------------------------------------------------------------------------------

# Establish location of home directory
docs = '/work2/vallender/Projects/Hall-Project/DocsAndFiles/'  # Set variable
a = docs  # Location of taxids.csv
os.chdir(docs)

#------------------------------------------------------------------------------
# Open taxids.csv which is a comma delimited list of all tax id's to be used.
#------------------------------------------------------------------------------

tax_ids = pd.read_csv('taxids.csv', header=None, dtype=str)  # 1st column = tax id's
tax_ids = list(tax_ids[0])


#------------------------------------------------------------------------------
# Create 'gi_lists' directory
#------------------------------------------------------------------------------

# Create directory for Gi lists
GiLists = a + "./gi_lists"  # Establish location of gi_lists Directory
b = GiLists  # Set a variable for the gi_lists directory
os.makedirs('%s' % b, exist_ok=True)  # Create the directory

#------------------------------------------------------------------------------
# Define the function we want to use.
#------------------------------------------------------------------------------

def get_gilists(ID):
    """ This function uses the blastdbcmd tool to get gi lists. It then uses the
    blastdb_aliastool to turn the list into a binary file."""
    # Use the accession #'s and the blastdbcmd tool to generate gi lists based on Organisms/Taxonomy ID's.
    os.system("blastdbcmd -db refseq_rna -entry all -outfmt '%g %T' | awk ' { if ($2 == " + ID + ") { print $1 } } ' > " + ID + "gi.txt")
    print(ID + "gi.txt has been created.")

    # Convert the .txt file to a binary file using the blastdb_aliastool.
    os.system("blastdb_aliastool -gi_file_in " + ID + "gi.txt -gi_file_out " + ID + "gi")
    print(ID + "gi binary file has been created.")

    # Remove the gi.text file
    os.system("rm -r *gi.txt")

    # Move the gi file to a folder
    os.system("mv " + ID + "gi gi_lists/")

#------------------------------------------------------------------------------
# if statement that gets to worker pool going to start multiprocessing
#------------------------------------------------------------------------------
procs = []
for tax in tax_ids:
    ID = tax
    proc = multiprocessing.Process(target=get_gilists, args=[ID])
    proc.start()
    procs.append(proc)

for proc in procs:
    proc.join()

print("Complete.")
