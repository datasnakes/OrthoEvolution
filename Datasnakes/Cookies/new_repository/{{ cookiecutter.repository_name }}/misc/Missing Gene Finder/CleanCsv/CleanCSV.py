# -*- coding: utf-8 -*-
"""
Last updated on February 15, 2017

@author: Shaurita D. Hutchins


This script is designed to remove duplicates from a .csv file and
not which duplicates were removed in a .txt file.
"""

# List of modules used.
import pandas as pd
import mygene as mg

#------------------------------------------------------------------------------
# Use pandas to read in the dataframe and create lists.

# Read in main file
maf = pd.read_csv('MAFV3.1.csv', index_col=False, dtype=str)

# Read in organisms file and create organisms list
orgs = pd.read_csv('Organisms.csv', index_col=False, dtype=str, header=None)
orglist = list(orgs[0])

#------------------------------------------------------------------------------
# Create dictionary for duplicate values and blank cells/values
dupdict = {}  # Dictionary of duplicates based on organisms/columns
nadict = {}  # Dictionary of black or n/a cells in the main file

# Create for loop that creates dicts
for org in orglist:
    dups = maf.duplicated(org, keep=False)
    dupdict[org] = dups
    nas = maf[org].isnull()
    nadict[org] = nas

#------------------------------------------------------------------------------
# Short definition that turns a dictionary into a csv file
def frametocsv(csvname, data):
    frame = pd.DataFrame.from_dict(data, orient='columns')
    frames = [maf.Tier, maf.Gene, frame]
    file = pd.concat(frames, axis=1)
    file.to_csv(csvname, index=False)

# Create the csv files
frametocsv(csvname='maf_duplicates_by_org.csv', data=dupdict)
frametocsv(csvname='maf_blanks.csv', data=nadict)