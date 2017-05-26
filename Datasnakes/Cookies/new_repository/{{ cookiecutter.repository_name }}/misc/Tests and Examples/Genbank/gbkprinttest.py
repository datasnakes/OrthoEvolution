# -*- coding: utf-8 -*-
"""
Created on Wed Sep 14 17:38:59 2016

@author: shutchins2

# Use one of the .gbk files in the Tests and Examples GitHub folder.
"""

from Bio import SeqIO
for index, record in enumerate(SeqIO.parse("HTR1A_Bos_taurus.gbk", "genbank")):
    print("index %i, ID = %s, length %i, with %i features"
          % (index, record.id, len(record.seq), len(record.features)) + "\n")
print(record + "\n")

dir(record)
