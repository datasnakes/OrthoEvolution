##############################################################################
# PyCharm Community Edition
# -*- coding: utf-8 -*-
"""
Orthologs-Project
missingORduplicate updated on 4/10/2017 at 11:27 AM
##############################################################################

    Input:

    Output:

    Description:

##############################################################################
@author: rgilmore
"""
##############################################################################
# Libraries:

import os
import pandas as pd
from os.path import dirname, abspath
from dir_mana import dir_mana
from lister import Lister

##############################################################################
# Custom Class Initializations
# :
# Use directory_management() class here so that we can stay organized
# and more easily access the proper directories on command
home = os.getcwd()
project = "Orthologs-Project"
user = "rgilmore"
where = dir_mana(home, project)
# Use lister() class here so that we can easily access our Master RNA Accession File
what = Lister('MAFV3.1.csv')  # Always make sure this file name is correct

## Add a path that contains custom libraries for import
#os.sys.path.append()
##############################################################################
# Global Initializations:

##############################################################################

# Read in main file
maf = pd.read_csv('data/processed/karg-maf.csv', index_col=False, dtype=str)

# Read in organisms file and create organisms list
orgs = pd.read_csv('data/interim/Organisms.csv', index_col=False, dtype=str, header=None)
orglist = list(orgs[0])

#------------------------------------------------------------------------------
# Create dictionary for duplicate values and blank cells/values
dupdict = {}  # Dictionary of duplicates based on organisms/columns
nadict = {}  # Dictionary of black or n/a cells in the main file

# Create for loop that creates dicts
for org in orglist:
    dups = maf.duplicated(org, keep=False)  # Get any duplicates in a column.
    dupdict[org] = dups
    nas = maf[org].isnull()  # Get any empty spaces in a column.
    nadict[org] = nas