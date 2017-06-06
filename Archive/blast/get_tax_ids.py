# -*- coding: utf-8 -*-
"""Get taxonomy ids from a list of organisms using Entrez."""
#------------------------------------------------------------------------------
# List of modules used.
#------------------------------------------------------------------------------

import os
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
docs = '/work2/vallender/Projects/Orthologs-Project/DocsAndFiles/'  # Set variable
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
else:
    print("Something is wrong.")

# Turn the list into a pandas dataframe and save it
df = pd.DataFrame(idlist, dtype=str)
df.to_csv('taxids.csv', header=False, index=False)
print('Your file taxids.csv has been saved.')

