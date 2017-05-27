# -*- coding: utf-8 -*-
"""
Created on Mon Sep 19 11:16:37 2016

@author: rgilmore
"""

import os
from Bio import SeqIO
import csv
from BioSQL import BioSeqDatabase
from fnmatch import fnmatch

server = BioSeqDatabase.open_database(driver="sqlite3", db="/work5/r2294/bin/NCBI_data/vertebrate_mammalian/DB/GPCR_Orthologs_DB.db")
count = 0
home = os.getcwd()
os.chdir('/work5/r2294/bin/NCBI_data/Raw_GBK_Files/HTR1A')
for file in (os.listdir('/work5/r2294/bin/NCBI_data/Raw_GBK_Files/HTR1A')):
    if fnmatch(file, '*.gbk') == False:
        continue
    
    print(file)
    db = server['HTR1A']
    try:
        c = db.load(SeqIO.parse(file, 'genbank'))
    except:
        print('error')
        continue
    server.commit()
    count = count +c
    print(count)
    input('db loaded')
os.chdir(home)
#count = 0
#for db_name in server.keys():
    #db = server[db_name]
#    for record in db.values():
#        for feature in record.features:
#            if feature.type == "source":
#                o = feature.qualifiers['organism']
#                print(feature.qualifiers['organism'][0])
#                print(o[0])
#            if feature.type == "CDS":
#                c = 0
#                #o = feature.qualifiers['organism']
#                #print(o[0])
#                input('there')
#                for x in feature.qualifiers['db_xref']:
#                    c += 1 
#                    if 'GI' in x:
#                        print(x)
#                        head, sup, tail = x.partition(':')
#                        print(head.lower())
#                        print('[' + tail)
#                        input('ok')
#                
#                input('ok')
#    count += 1
#    print(count)
db = server['HTR2A']
for record in db.values():
    db_name = 'HTR1A'
    print(record.format('genbank'))
    print(dir(record.annotations))
    print(record.annotations)
    for x in record.annotations:
        print(x.values)
print(count, 'Number of files')
try:
    record = db.lookup(accession = 'XM_007486379')
    print(record.format('genbank'))
    print(record.name, 'name')
    print(record.id, 'id')
    print(record)
except IndexError:
    print('XM_007486379.2 not in %s' % db_name)
#    continue
