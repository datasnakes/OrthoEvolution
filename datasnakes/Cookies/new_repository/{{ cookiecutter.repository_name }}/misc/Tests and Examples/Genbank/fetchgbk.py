# -*- coding: utf-8 -*-
"""
Date created: Wed Mar  1 10:39:46 2017
Author: S. Hutchins

Script description: Quick script to retrieve genbank files and extract features

"""

# Modules used
from Bio import Entrez
import time as t
import pandas as pd
from Bio import SeqIO
import os

# Change to gbk dir
gbkdir = '/work2/vallender/Projects/Hall-Project/GbkDir'
os.chdir(gbkdir)

# Always tell NCBI who you are
Entrez.email = "shutchins2@umc.edu"

# Create a short variable for esearch, Entrez.read, and efetch
efetch = Entrez.efetch

# Import the csv file
file = pd.read_csv('/work2/vallender/Projects/Hall-Project/DocsAndFiles/Master_GI_File.csv', header=None)

# Transpose the dataframe
df = pd.DataFrame(file, index=None)
df = df.T
df = df.fillna("none")
df = df.ix[1:]  # Remove the first row

idlist = list(df[1])
orglist = list(df[0])

# Retrieve a genbank file using the GI number
for ID, org in zip(idlist, orglist):
    if ID == "none":
        pass
    else:
        nucl = efetch(db="nucleotide", id=ID, rettype="gb", retmode="text")
        gbkfile = open("APOL1_" + org + ".gbk", "w")
        gbkfile.write(nucl.read())
        nucl.close()
        gbkfile.close()
        record = SeqIO.read("APOL1_" + org + ".gbk", "genbank")

        # Change to cds dir
        cdsdir = '/work2/vallender/Projects/Hall-Project/CdsDir'
        os.chdir(cdsdir)
        output_handle = open("APOL1_" + org + "_cds_nucl.fasta", "w")

        # For loop to extract features
        for feature in record.features:
                # Other annotated features are 'Gene', 'mRNA', 'CDS', and 'ncRNA'.
                if feature.type == "CDS":
                    # Use record.dbxrefs here. Look up record features in Ipython
                    # using 'dir(record)'.
                    feature_name = org
                    feature_seq = feature.extract(record.seq)
                    # Simple FASTA output without line wrapping:
                    output_handle.write(
                        ">" + feature_name + "  " + "\n" + str(feature_seq) + "\n")
                    output_handle.close()
        t.sleep(10)