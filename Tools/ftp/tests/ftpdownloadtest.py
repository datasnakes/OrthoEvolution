# -*- coding: utf-8 -*-
"""
Date created: Tue Feb 28 17:32:04 2017
Author: S. Hutchins

Script description:

"""

# Modules used
from ftplib import FTP, error_perm
import os
import fnmatch

#------------------------------------------------------------------------------
# Connect to the NCBI directory
ncbi = 'ftp.ncbi.nlm.nih.gov/'
refseqrna = '/blast/db/'  # Connect to this blastdb
blast = '/blast/'
ftp = FTP("ftp.ncbi.nlm.nih.gov", timeout=None)

# Login using email as password
ftp.login(user='anonymous', passwd='shutchins2@umc.edu')

# Change to the desired directory
#ftp.cwd(refseqrna)
ftp.cwd(blast)
# Use ftp.pwd() to find out the current directory

# Get a list of all the files in the directory
#ftp.retrlines('LIST')

# This is a list of the file names
filenames = ftp.nlst()

# Create a for loop that downloads the files
for filename in filenames:
    if fnmatch.fnmatch(filename, 'blastftp*'):
        host_file = os.path.join(filename)
        try:
            with open(host_file, 'wb') as local_file:
                ftp.retrbinary('RETR %s' % filename, local_file.write)
        except error_perm:
            pass

ftp.quit()