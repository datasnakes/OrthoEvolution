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
Project Name: Orthologs Project
"""
# Modules used
from ftplib import FTP, error_perm
import os
import fnmatch
import sys
import logging as log
from datetime import datetime as d
import subprocess

#------------------------------------------------------------------------------
# Set up the logger for logging
format1 = '%a %b %d %I:%M:%S %p %Y'  # Used to add as a date
format2 = '%m-%d-%Y@%I:%M:%S-%p'  # Used to append to archives
format3 = '%m-%d-%Y'

log.basicConfig(filename="logs/blast_setup_" + d.now().strftime(format3) + ".log", level=log.INFO)
log.info("#------------------------------------------------------------------")
log.info("The script name is %s" % os.path.basename(__file__))
log.info("The date and time is currently %s" % str(d.now().strftime(format1)))
log.info("#------------------------------------------------------------------")
log.info("Download and update the refseq_rna blast database.")

#------------------------------------------------------------------------------
home = os.getcwd()  # My home directory for this script/project
dbpath = '/work5/r2295/bin/databases/refseq_rna_db'  # My current dbpath

# Create a directory for the database if one doesn't exist
try:
    # If the directory exists,
    if os.path.exists(dbpath) == True:
        log.info("The %s directory exists & will be archived." % str(dbpath))
        # Move any files that are in the directory to a dated archive folder.
        # Moving a directory in linux/unix essentially renames it.
        os.system('mv ' + dbpath + ' /work5/r2295/bin/databases/refseqrnadb_archive_' + d.now().strftime(format2))
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
    else:
        log.info("Something went wrong.")
        ftp.quit()
except os.error:
    log.info("There has been an error.")
    log.info(sys.exit())
    log.info(ftp.quit())

ftp.close()

#------------------------------------------------------------------------------
# Unzip all of the files and remove unneccessary files
os.system("for file in *.tar.gz; do tar xvf $file; done")  # Unzip the database files
os.system("rm -r *.tar.gz")
log.info("The files have been unzipped, and Part 1 has finished.")
log.info("#------------------------------------------------------------------")
