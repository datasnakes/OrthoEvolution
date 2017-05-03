# -*- coding: utf-8 -*-
"""
File Name: blast_setup.py
Description: This script 1) updates/downloads the refseq_rna blast db files,
2) creates a list of taxonomy ids based on the list of organisms, and 3) create
a csv file with only human accessions and genes for downstream usage. Check the
ReadMe file for a first time setup of a blast database. Both steps occur on the
head node of the MCSR because they need internet access.

@author: Shaurita Hutchins
Date Created: Tue Feb 28 19:05:41 2017
Project Name: KARG Project
"""
# Modules used
from ftplib import FTP, error_perm
import os
import fnmatch
import sys
import logging as log
from datetime import datetime as d
from Bio import Entrez
from time import sleep as pause
import pandas as pd
import subprocess

#------------------------------------------------------------------------------
# Set up the logger for logging
format1 = '%a %b %d %I:%M:%S %p %Y'  # Used to add as a date
format2 = '%m-%d-%Y@%I:%M:%S-%p'  # Used to append to archives
format3 = '%m-%d-%Y'  # Use for Windows

log.basicConfig(filename="blast_setup_" + d.now().strftime(format3) + ".log", level=log.INFO)
log.info("#------------------------------------------------------------------")
log.info("The script name is %s" % os.path.basename(__file__))
log.info("The date and time is currently %s" % str(d.now().strftime(format1)))
log.info("#------------------------------------------------------------------")
log.info("PART 1")
log.info("Download and update the refseq_rna blast database.")

#------------------------------------------------------------------------------
home = os.getcwd()  # My home directory for this script/project
dbpath = 'databases/refseq_rna_db'  # My current dbpath

# Create a directory for the database if one doesn't exist
try:
    # If the directory exists,
    if os.path.exists(dbpath) == True:
        log.info("The %s directory exists & will be archived." % str(dbpath))
        # Move any files that are in the directory to a dated archive folder.
        # Moving a directory in linux/unix essentially renames it.
        os.system('mv ' + dbpath + ' /databases/refseqrnadb_archive_' + d.now().strftime(format2))
        os.mkdir(dbpath)   # Recreate the database directory
        log.info("The %s directory was created." % str(dbpath))
        pass
    else:  # If the directory does not exist
        os.mkdir(dbpath)
        log.info("The %s directory was created." % str(dbpath))
except os.error:
    log.info("Error.")
    sys.exit("There has been an error.")

os.chdir(dbpath)  # Change to the database directory

#------------------------------------------------------------------------------
# Connect to the NCBI ftp site
try:
    ncbi = 'ftp.ncbi.nlm.nih.gov/'
    blastdb = '/blast/db/'  # Set variable for the blastdb subdirectory
    ftp = FTP("ftp.ncbi.nlm.nih.gov", timeout=None)
    # Login using email as password
    ftp.login(user='anonymous', passwd='shutchins2@umc.edu')
    log.info("Successful FTP login.")
except error_perm:  # This error will be thrown if there's a connection issue.
    log.info("FTP connection error.")
    sys.exit()

# Change to the desired directory
ftp.cwd(blastdb)
# Use ftp.pwd() to find out the current directory
# Use ftp.retrlines('LIST') to get a list of all the files in the directory

# This is a list of the file names in the current directory
filenames = ftp.nlst()

#------------------------------------------------------------------------------
# Create a for loop that writes the list/text file of files wanted
with open('downloadlist.txt', 'w') as downloads:
    for filename in filenames:
        if fnmatch.fnmatch(filename, 'refseq_rna*'):  # Get only those files.
            refseq_file = os.path.join(filename)
            # Write the url of each refseq_rna db file to a text file.
            downloads.writelines(ncbi + blastdb + refseq_file + '\n')
        # use elif here to get the taxdb.tar.gz file.
        elif fnmatch.fnmatch(filename, 'taxdb*'):
            taxdb_file = os.path.join(filename)
            downloads.writelines(ncbi + blastdb + taxdb_file + '\n')

# Download the list of files using 'wget' on linux/unix
try:
    cmd = 'cat downloadlist.txt | xargs -n 1 -P 8 wget'
    status = subprocess.call([cmd], shell=True)
    if status == 0:
        log.info("The refseqrna blast db files have downloaded.")
        # Unzip all of the files and remove unneccessary files
        os.system("for file in *.tar.gz; do tar xvf $file; done")  # Unzip the database files
        os.system("rm -r *.tar.gz")
        log.info("The files have been unzipped, and Part 1 has finished.")
        log.info("#------------------------------------------------------------------")
    else:
        log.info("Something went wrong.")
        ftp.quit()
except os.error:
    log.info("There has been an error.")
    log.info(sys.exit())
    log.info(ftp.quit())

ftp.close()
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------

# Mark start of PART 2 with logging information
log.info("#------------------------------------------------------------------")
log.info("PART 2")
log.info("Create a list of taxonomy ids based on the organisms list.")

#------------------------------------------------------------------------------
# Change to the home directory
os.chdir(home)

#------------------------------------------------------------------------------
# Open Orgnanisms.csv which is a comma delimited list of all species names to be used.
# Read .csv into pandas. If no columns for the data, set header as None or 0 (int)
# If the data has columns, simply just read in the csv without the header option
orgsfile = pd.read_csv('data/initial-data/organisms.csv', header=None)
# Create a list name/variable and use list()
orgs = list(orgsfile[0])

#------------------------------------------------------------------------------
# Get a list of taxonomy ids for a corresponding organisms list
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

log.info("This is the id list:  %s" % idlist)
if len(idlist) == len(orgs):
    log.info('There are the same number of taxonomy ids as organisms.')
else:
    log.info("Something is wrong.")

# Turn the list into a pandas dataframe and save it
df = pd.DataFrame(idlist, dtype=str)
df.to_csv('data/initial-data/taxids.csv', header=False, index=False)
log.info('Your file taxids.csv has been saved.')

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------

# Mark start of PART 3 with logging information
log.info("#------------------------------------------------------------------")
log.info("PART 3")
log.info("Create a csv file with only human accessions and genes for downstream usage.")

#------------------------------------------------------------------------------
# Open the file with the list of all the genes, accessions, and tiers
mal = pd.read_csv('data/initial-data/master_accessions_list.csv')
human_acc = list(mal.Homo_sapiens)  # Create a list of human accessions from column
genes_list = list(mal.Gene)  # Create a list of genes from the column
tiers_list = list(mal.Tier)

# Turn the lists into dataframes
accs = pd.DataFrame(human_acc, dtype=str)
genes = pd.DataFrame(genes_list, dtype=str)
tiers = pd.DataFrame(tiers_list, dtype=str)

# List of dataframes I want to combine
frames = [tiers, genes, accs]

# Merge the dataframes into 1 dataframe
human_acc_info = pd.concat(frames, axis=1)

# Save the dataframe as a csv file which will be used in AccCollect.py
human_acc_info.to_csv('data/initial-data/homo_sapiens_accessions.csv', index=False, header=None)

log.info("homo_sapiens_accessions.csv file was created.")
log.info("The script has finished.")
log.info("#------------------------------------------------------------------")
log.info("#------------------------------------------------------------------")
log.info("#------------------------------------------------------------------")

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------