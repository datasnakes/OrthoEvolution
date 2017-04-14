##############################################################################
# PyCharm Community Edition
# -*- coding: utf-8 -*-
"""
GPCR-Orthologs-Project
dataviz updated on 1/10/2017 at 6:27 PM
##############################################################################

    Input:

    Output:

    Description:

##############################################################################
@author: Work
"""
##############################################################################
# Libraries:

import os

from BioSQL import BioSeqDatabase

from dir_mana import dir_mana
from lister import Lister

##############################################################################
# Custom Class Initializations
# :
# Use directory_management() class here so that we can stay organized
# and more easily access the proper directories on command
home = os.getcwd()
project = "GPCR-Orthologs-Project"
user = "Work"
where = dir_mana(home, project)
# Use lister() class here so that we can easily access our Master RNA Accession File
what = Lister('MAFV3.1.csv')  # Always make sure this file name is correct

## Add a path that contains custom libraries for import
#os.sys.path.append()
##############################################################################
# Global Initializations:

##############################################################################

# genes = '%s' % genes.r_repr().rstrip('\n')
arg_string1 = 'gene'
db_list = []
for gene in what.master_gene_list:
    arg_string1 += (' %s' % gene)
arg_string2 = 'org'
for org in what.org_list:
    arg_string2 += (' %s' % org)
for file in os.listdir('/srv/shiny-server/GPCR-Orthologs-Project/CODE/1_Databases/Vallender Data/'):
    if '.db' in str(file):
        db_list.append(file)
for db_name in db_list:
    server = BioSeqDatabase.open_database(driver='sqlite3', db='/srv/shiny-server/GPCR-Orthologs-Project/CODE/1_Databases/Vallender Data/%s' % db_name)
    for sub_db_name in server.keys():
        db = server[sub_db_name]
        Accession = str(what.A_Get1['HTR1A']['Homo_sapiens'])
        Accession, Sup, Version = Accession.partition('.')
        Accession = Accession.upper()
        try:
            record = db.lookup(accession=Accession)
            with open('/srv/shiny-server/GPCR-Orthologs-Project/CODE/1_Databases/Vallender Data/HTR1A_Homo_sapiens.gbk', 'w') as GB_file:
                GB_file.write(record.format('genbank'))
                print(GB_file.name, 'created')
            server_flag = True
        except IndexError:
            continue

os.system('Rscript /srv/shiny-server/GPCR-Orthologs-Project/CODE/Lib/web.R %s %s'
          % (arg_string1, arg_string2))

#//TODO-ROB Add code that takes the features and turns them into actual sequence data
#//TODO-ROB Add code that breaks up the annotation attribute
