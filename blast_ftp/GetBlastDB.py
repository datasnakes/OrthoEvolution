# -*- coding: utf-8 -*-
"""
Date created: Tue Feb 28 19:05:41 2017
Author: S. Hutchins

Script description: This script downloads the refseq_rna blast db files.

"""

# Modules used
from ftplib import FTP, error_perm
import os
import fnmatch
# import logging as log

#------------------------------------------------------------------------------
# Connect to the NCBI directory
ncbi = 'ftp.ncbi.nlm.nih.gov/'
refseqrna = '/blast/db/'  # Connect to this blastdb
ftp = FTP("ftp.ncbi.nlm.nih.gov", timeout=None)

# Login using email as password
ftp.login(user='anonymous', passwd='shutchins2@umc.edu')

# Change to the desired directory
ftp.cwd(refseqrna)
# Use ftp.pwd() to find out the current directory

# Get a list of all the files in the directory
ftp.retrlines('LIST')

# This is a list of the file names
filenames = ftp.nlst()

# Create a for loop that downloads the files
for filename in filenames:
    if fnmatch.fnmatch(filename, 'refseq_rna*'):
        host_file = os.path.join(filename)
        try:
            with open(host_file, 'wb') as local_file:
                ftp.retrbinary('RETR %s' % filename, local_file.write)
                print("The following file has been downloaded: " + filename)
        except error_perm:
            print("There has been an error.")
            pass

ftp.quit()