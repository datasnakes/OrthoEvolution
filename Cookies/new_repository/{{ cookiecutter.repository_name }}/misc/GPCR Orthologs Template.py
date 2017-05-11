#set($hashtags = '##############################################################################')
${hashtags}
# ${PRODUCT_NAME}
# -*- coding: utf-8 -*-
"""
${PROJECT_NAME}
${NAME} updated on ${DATE} at ${TIME}
${hashtags}

    Input:

    Output:

    Description:

${hashtags}
@author: ${USER}
"""
${hashtags}
# Libraries:


import os
from os.path import dirname, abspath
${hashtags}
# Directory Initializations:

# All of our static directories are initialized here so that we can stay organized
# and more easily access the proper directories on command

CWD = os.getcwd() # */GPCR-Orthologs-Project/.idea/fileTemplates
_homedir = dirname(dirname(abspath(__file__))) # */GPCR-Orthologs-Project
_CODEdir = _homedir + '/CODE' # */GPCR-Orthologs-Project/CODE
_DBdir = _CODEdir + '/1 Databases' # */GPCR-Orthologs-Project/CODE/1 Databases
_AGdir = _CODEdir + '/2_Accession_Gathering' # */GPCR-Orthologs-Project/CODE/2_Accession_Gathering
_GB_FASTAdir = _CODEdir + '/3_GenBank-FASTA' # */GPCR-Orthologs-Project/CODE/3_GenBank-FASTA
_ALIGNdir = _CODEdir + '/4 Alignment' # */GPCR-Orthologs-Project/CODE/4 Alignment
_PHYLOdir = _CODEdir + '/5 Phylogenetic Analysis' # */GPCR-Orthologs-Project/CODE/5 Phylogenetic Analysis
_MISCdir = _CODEdir + '/Miscellaneous' # */GPCR-Orthologs-Project/CODE/Miscellaneous

#set($c_Lib = '## CUSTOM LIBRARIES')
${c_Lib}
_LIBdir = _CODEdir  + '/Lib' # */GPCR-Orthologs-Project/CODE/Lib


#set($r1 = '## Add a path that contains custom libraries for import')
${r1}
os.sys.path.append(_LIBdir)
${hashtags}
# Initializations:

${hashtags}