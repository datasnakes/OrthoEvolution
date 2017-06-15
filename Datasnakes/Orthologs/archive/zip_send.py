# -*- coding: utf-8 -*-
"""
File Name:
Description:

Author: shutchins2
Date Created: Tue Apr 11 12:31:17 2017
Project Name:
"""

import os
import logging as log
import pandas as pd
from datetime import datetime as d
import zipfile
import pexpect

#------------------------------------------------------------------------------
# Set up logging
log.basicConfig(filename="logs/zip_send.log", level=log.INFO)
log.info("#------------------------------------------------------------------")
log.info("The script name is %s" % os.path.basename(__file__))
log.info("The date and time is currently %s" % str(d.now()))
log.info("#------------------------------------------------------------------")
log.info("Run clustal omega.")

#------------------------------------------------------------------------------
# Create the main output directories
# Home Directory
home = os.getcwd() + '/'
h = home

clustal_out = 'data/clustal-output/'
paml_out = 'data/paml-output/'
phyml_out = 'data/phyml-output/'
processed = 'data/processed/'

dir_list = [clustal_out, phyml_out, paml_out, processed]  # List of directories

for directory in dir_list:
    if os.path.exists(directory) == True:
        log.info('The directory %s exists.' % directory)
    else:
        os.mkdir(directory)
        log.info('The directory %s has been created' % directory)

#------------------------------------------------------------------------------
g = pd.read_csv(clustal_out + 'genes_to_align.txt', sep='\t', header=None)
geneslist = list(g[0])

#------------------------------------------------------------------------------
fileslist = []
for gene in geneslist:
    os.chdir(clustal_out)
    clustfile = gene + "_aligned_cds_nucl.fasta"
    clustpath = gene + "_Aligned/" + clustfile
    if os.path.exists(clustpath) == True:
        log.info("%s exists. Append it to the files list." % clustfile)
        fileslist.append(clustfile)
        os.system('cp ' + clustpath + ' ' + h + processed + 'files-to-zip/')
        os.chdir(h)
    else:
        log.info("%s does not exist." % clustfile)
        os.chdir(h)
        pass
#------------------------------------------------------------------------------
# Save the list of files that will be zipped
df = pd.DataFrame(fileslist)
df.to_csv(
    processed +
    'files-to-zip/alignment_files_to_zip.txt',
    sep='\t',
    index=False,
    header=None)

#------------------------------------------------------------------------------
# Send all alignment files to the server.
os.chdir(processed + 'files-to-zip')
sendfile = pexpect.spawnu(
    "scp *.fasta shutchins2@162.243.56.106:/srv/shiny-server/public-ftp/data/karg/alignments")
sendfile.sendline("shutchins2\r")
sendfile.waitnoecho()

#------------------------------------------------------------------------------
# Create and save the zip file
with zipfile.ZipFile(h + processed + 'karg_alignments.zip', 'w') as clustzip:
    for f in fileslist:
        clustzip.write(f)
        log.info("%s has been written to the zip file." % f)

log.info("The zip file (%s) has been created and saved." % clustzip)

#------------------------------------------------------------------------------
