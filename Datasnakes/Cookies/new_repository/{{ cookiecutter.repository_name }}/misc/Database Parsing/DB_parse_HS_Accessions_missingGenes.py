# -*- coding: utf-8 -*-
"""
Created on Fri Sep 23 16:07:00 2016

@author: rgilmore
"""
import sqlite3
import fnmatch
import time
import Bio
import csv
import os
import pexpect
from Bio import SeqIO
from BioSQL import BioSeqDatabase


server_db_list = [] #list of db names.  I will create a DB for the orthologs of interest soon.
Acc0_list = []
Acc1_list = []
dbb = []
s_count = 0
o_count = 0
#Open the server that we want to look at
server = BioSeqDatabase.open_database(driver="sqlite3", db="/work5/r2294/bin/NCBI_data/vertebrate_mammalian/DB/GPCR_Orthologs_DB.db")
#Get the current working directory and set it to the home variable
home = os.getcwd()
with open('Org_list.csv', 'r') as csvfile:
    O = csv.reader(csvfile, delimiter=',')
    org_list = []
    for line in O:
        for cell in line:
            o_count += 1
            org_list.append(cell)
    print(o_count)
#make a list of database names
for db_name in server.keys():
    s_count += 1
    print(s_count)
    server_db_list.append(db_name)
    if db_name in org_list:
        print(db_name)
print(server_db_list)
input('this is the server_db_list Do you enjoy it?  Very much eh? ......')
with open('DB_HS_Accessions.csv', 'r') as csvfile:
    Acc = csv.reader(csvfile)
    for line in Acc:
        c_count = 0
        for cell in line:
            c_count +=1
            if c_count == 1 or c_count == 3:
                continue
            if cell not in org_list:
                print(cell)
#To start you have to parse each subdatabase on the "server" in order to search through each one.
for db_name in server.keys():
    dbb.append(db_name)
    org_flag0 = False
    db = server[db_name]
    print(db_name)
    for record in db.values():
        if org_flag0 == True:
            break
        for feature in record.features:
            if feature.type == 'source':
                print(feature.qualifiers['organism'][0])
                if feature.qualifiers['organism'][0] == 'Homo sapiens':
                    Acc0_list.append(feature.qualifiers['organism'][0])
                    org_flag0 = True
                    print('\n\n\n', 'Here it is:')
                    print(feature.qualifiers['organism'][0])
                    print(record.id)
                    
                    print('qualifiers count: ', len(Acc0_list))
                    
                    if org_flag0 == True:
                        break
        if 'Homo sapiens' in record.description:
            Acc1_list.append(db_name)
            print('description count: ', len(Acc1_list))
    print('db count: ', len(dbb))
    if len(Acc1_list) != len(dbb):
        input('stop')
        dbb.pop(0)
        
        
            