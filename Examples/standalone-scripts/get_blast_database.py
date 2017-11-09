"""This standalone script downloads a NCBI ftp databases."""
from ftplib import FTP, error_perm
import os
import fnmatch
import sys
import logging as log
from datetime import datetime as d
import subprocess
import contextlib

# TODO import our ftp client.
# TODO import our logging

# Set up the logger for logging
format1 = '%a %b %d %I:%M:%S %p %Y'  # Used to add as a date
format2 = '%m-%d-%Y@%I:%M:%S-%p'  # Used to append to archives
format3 = '%m-%d-%Y'

home = os.getcwd()  # My home directory for this script/project
dbpath = home  # My current dbpath

# Create a directory for the database if one doesn't exist
try:
    # If the directory exists,
    if os.path.exists(dbpath):
        # Move any files that are in the directory to a dated archive folder.
        # Moving a directory in linux/unix essentially renames it.
        os.system(
            'mv ' +
            dbpath +
            ' /work5/r2295/bin/databases/refseqrnadb_archive_' +
            d.now().strftime(format2))
        os.mkdir(dbpath)   # Recreate the database directory
        pass
    else:  # If the directory does not exist
        os.mkdir(dbpath)
        log.info("The %s directory was created." % str(dbpath))
except os.error:
    log.info("Error.")
    sys.exit("There has been an error.")

os.chdir(dbpath)  # Change to the database directory

# Connect to the NCBI ftp site
try:
    ncbi = 'ftp.ncbi.nlm.nih.gov/'
    blastdb_path = '/blast/db/'  # Set variable for the blastdb subdirectory
    ftp = FTP("ftp.ncbi.nlm.nih.gov", timeout=None)
    # Login using email as password
    ftp.login(user='anonymous', passwd='shutchins2@umc.edu')
    log.info("Successful FTP login.")
except error_perm:  # This error will be thrown if there's a connection issue.
    log.info("FTP connection error.")
    sys.exit()

# Change to the desired directory
ftp.cwd(blastdb_path)
# Use ftp.pwd() to find out the current directory
# Use ftp.retrlines('LIST') to get a list of all the files in the directory

# This is a list of the file names in the current directory
filenames = ftp.nlst()

# Create a for loop that writes the list/text file of files wanted
with open('downloadlist.txt', 'w') as downloads:
    for filename in filenames:
        if fnmatch.fnmatch(filename, 'refseq_rna*'):  # Get only those files.
            refseq_file = os.path.join(filename)
            # Write the url of each refseq_rna db file to a text file.
            downloads.writelines(ncbi + blastdb_path + refseq_file + '\n')
        # use elif here to get the taxdb.tar.gz file.
        elif fnmatch.fnmatch(filename, 'taxdb*'):
            taxdb_file = os.path.join(filename)
            downloads.writelines(ncbi + blastdb_path + taxdb_file + '\n')

# Download the list of files using 'wget' on linux/unix
with contextlib.suppress(os.error):
    cmd = 'cat downloadlist.txt | xargs -n 1 -P 8 wget'
    status = subprocess.call([cmd], shell=True)
    if status == 0:
        log.info("The refseqrna blast db files have downloaded.")
    else:
        log.info("Something went wrong.")
        ftp.quit()

ftp.close()

# Unzip all of the files and remove unneccessary files
# Unzip the database files
os.system("for file in *.tar.gz; do tar xvf $file; done")
os.system("rm -r *.tar.gz")
log.info("The files have been unzipped, and Part 1 has finished.")
log.info("#------------------------------------------------------------------")