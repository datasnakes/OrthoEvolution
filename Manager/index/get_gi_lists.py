# -*- coding: utf-8 -*-
"""
File Name: get_gi_lists.py
Description: This script is designed to create a gi list based on the refseq_rna database
for each taxonomy id on the MCSR. It will also convert the gi list into a
binary file which is more efficient to use with NCBI's Standalone Blast tools.

@author: Shaurita D. Hutchins
Date Created: Wed Mar 29 16:40:01 2017
Project Name: Addictions Project
"""
# Modules used
import os
from datetime import datetime as d
import pandas as pd
from time import time
from multiprocessing import Pool
import logging as log

# Set up logging
log.basicConfig(filename="get_gi_lists.log", level=log.INFO)
log.info("#------------------------------------------------------------------")
log.info("The script name is %s" % os.path.basename(__file__))
log.info("The date and time is currently %s" % str(d.now()))
log.info("#------------------------------------------------------------------")
log.info("Create a GI list for each organism using the taxonomy id and the blastdbcmd tool on the MCSR.")

# Get a taxid list from the file and change to the output directory
tax_ids = pd.read_csv('taxids.csv', header=None, dtype=str)  # 1st column = tax id's
tax_ids = list(tax_ids[0])
os.chdir('data')


def get_gilists(id):
    """ This function uses the blastdbcmd tool to get gi lists. It then uses the
    blastdb_aliastool to turn the list into a binary file.
    The input (id) for the function is a taxonomy id.
    """
    # Use the accession #'s and the blastdbcmd tool to generate gi lists based on Organisms/Taxonomy id's.
    os.system("blastdbcmd -db refseq_rna -entry all -outfmt '%g %T' | awk ' { if ($2 == " + id + ") { print $1 } } ' > " + id + "gi.txt")
    log.info(id + "gi.txt has been created.")
    # Convert the .txt file to a binary file using the blastdb_aliastool.
    os.system("blastdb_aliastool -gi_file_in " + id + "gi.txt -gi_file_out " + id + "gi")
    log.info(id + "gi binary file has been created.")
    # Remove the gi.text file
    os.system("rm " + id + "gi.txt")
    log.info(id + "gi.text file has been deleted.")


def main(idlist):
    """This function uses a pool to start multiple processes to get gi lists.
    The argument (idlist) should be a list of taxonomy ids or 1 id.
    """
    ts = time()
    with Pool(processes=10) as p:
        p.map(get_gilists, idlist)
        log.info("Took {} minutes to get all gi lists.".format((time() - ts)/60))
# ------------------------------------------------------------------------------
# Run the main function that creates the lists in parallel.
main(idlist=tax_ids)
