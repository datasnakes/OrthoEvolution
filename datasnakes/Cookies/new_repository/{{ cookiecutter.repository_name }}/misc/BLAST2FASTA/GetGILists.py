# -*- coding: utf-8 -*-
"""
Last updated on January 27, 2017

@author: Shaurita D. Hutchins


This script is designed to create a gi list based on the refseq_rna database
for each taxonomy id on the MCSR. It will also convert the gi list into a
binary file which is more efficient to use with NCBI's Standalone Blast tools.

"""

#------------------------------------------------------------------------------
# List of modules used.
#------------------------------------------------------------------------------

import os
import csv
import sys
from Bio import Entrez
from time import sleep as pause
import pandas as pd


#------------------------------------------------------------------------------
# Mark start of program with printed text.
#------------------------------------------------------------------------------

print("\n" + (90 * "#") + "\n" + "Create a GI list for each organism using the \
taxonomy id and the blastdbcmd tool on the MCSR.  ####" + "\n" + (90 * "#") + "\n")

#------------------------------------------------------------------------------
# Establish home directory
#------------------------------------------------------------------------------

# Establish location of home directory
docs = '/work2/vallender/Projects/GPCR-Orthologs/DocsAndFiles/'  # Set variable
a = docs  # Location of taxids.csv
os.chdir(docs)
#------------------------------------------------------------------------------
# Open Orgnanisms.csv which is a comma delimited list of all species names to be used.
#------------------------------------------------------------------------------

# Read .csv into pandas. If no columns for the data, set header as None or 0 (int)
# If the data has columns, simply just read in the csv without the header option
orgs = pd.read_csv('Organisms.csv', header=None)
# Create a list name/variable and use list()
orgs = list(orgs[0])
#------------------------------------------------------------------------------
# Get a list of taxonomy ids for a corresponding organisms list
#------------------------------------------------------------------------------
""" One thing I will add later is creating a list of organisms to make sure
that the length of both taxid list and organism list are the same.

Could also add - using pandas - the creation of a file with both the taxid list
and organism list.
"""

# Always tell ncbi whom you are
Entrez.email = "shutchins2@umc.edu"

# Create a short variable for esearch & record reading
esearch = Entrez.esearch
read = Entrez.read

# Initialize a list
idlist = []

# Create a for loop to compile a list of taxonomy ids
for species in orgs:

    # Create a handle for the search
    searchresult = esearch(db="taxonomy", term=species)

    # Read and print the record or results of the search
    record = read(searchresult)

    # Replace unwanted characters
    idresult = record['IdList']
    x = str(idresult)
    x = x.replace("'", "")
    x = x.replace("[", "")
    x = x.replace("]", "")

    # Append the id to the id list
    idlist.append(x)

    # Keep updated on how this is progressing
    print("Taxonomy ID: " + x, "| Species Name: " + species)
    # Make a short pause (2 seconds) as to not break the entrez rules
    pause(2)

print(idlist)
if len(idlist) == len(orgs):
    print('There are the same number of taxonomy ids as organisms.')

# Turn the list into a pandas dataframe and save it
df = pd.DataFrame(idlist, dtype=str)
df.to_csv('taxids.csv', header=False, index=False)
print('Your file taxids.csv has been saved.')
pause(10)

#------------------------------------------------------------------------------
# Open taxids.csv which is a comma delimited list of all tax id's to be used.
#------------------------------------------------------------------------------

taxid = open('taxids.csv')  # 1st column = tax id's
t = csv.reader(taxid)
TaxID_count = 0


#------------------------------------------------------------------------------
# Create 'gi_lists' directory
#------------------------------------------------------------------------------

# Create directory for Gi lists
GiLists = a + "./gi_lists"  # Establish location of gi_lists Directory
b = GiLists  # Set a variable for the gi_lists directory
os.makedirs('%s' % b, exist_ok=True)  # Create the directory

#------------------------------------------------------------------------------
# Main "for" loop in this script that creates gi list files in binary format
# in the current/home directory.
#------------------------------------------------------------------------------

for ID in t:
    TaxID_count = TaxID_count + 1

#------------------------------------------------------------------------------
    # Use the accession #'s and the blastdbcmd tool to generate gi lists based on Organisms/Taxonomy ID's.
    os.system("blastdbcmd -db refseq_rna -entry all -outfmt '%g %T' | awk ' { if ($2 == " + str(ID[0]) + ") { print $1 } } ' > " + str(ID[0]) + "gi.txt")
    print(str(ID[0]) + "gi.txt" + " has been created.")

#------------------------------------------------------------------------------
    # Convert the .txt file to a binary file using the blastdb_aliastool.
    os.system("blastdb_aliastool -gi_file_in " + str(ID[0]) + "gi.txt -gi_file_out " + str(ID[0]) + "gi")
    print(str(ID[0]) + "gi binary file has been created.")

#------------------------------------------------------------------------------
    # Remove the gi.text file
    os.system("rm -r *gi.txt")

#------------------------------------------------------------------------------
    # Move the gi file to a folder
    os.system("mv " + str(ID[0]) + "gi gi_lists/")

sys.exit("The script has completed! ✓")
