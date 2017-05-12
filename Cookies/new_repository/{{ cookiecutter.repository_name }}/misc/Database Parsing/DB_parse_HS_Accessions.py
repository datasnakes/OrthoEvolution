# -*- coding: utf-8 -*-
"""
Created on Wed Sep 14 12:27:33 2016

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

org_flag0 = False
org_flag1 = False
org_flag2 = False
org_flag3 = False

server_db_list = [] #list of db names.  I will create a DB for the orthologs of interest soon.
Acc_list = []
#Open the server that we want to look at
server = BioSeqDatabase.open_database(driver="sqlite3", db="/work5/r2294/bin/NCBI_data/vertebrate_mammalian/DB/GPCR_Orthologs_DB.db")
#Get the current working directory and set it to the home variable
home = os.getcwd()

#make a list of database names
for db_name in server.keys():
    server_db_list.append(db_name)
print(server_db_list)
input('this is the server_db_list Do you enjoy it?  Very much eh? ......')
with open('DB_HS_Accessions.csv', 'w') as csvfile:
    Acc = csv.writer(csvfile)
    Acc.writerow(['Organism', 'Gene', 'Accession'])
#To start you have to parse each subdatabase on the "server" in order to search through each one.
for db_name in server.keys():
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
                    Acc_list.append(feature.qualifiers['organism'][0])
                    org_flag0 = True
                    with open("DB_HS_Accessions.csv", 'a') as csvfile:
                        Acc_list.append(db_name)
                        Acc_list.append(record.id)
                        Acc = csv.writer(csvfile)
                        print('\n\n\n', 'Here it is:')
                        print(feature.qualifiers['organism'][0])
                        print(record.id)
                        Acc.writerow(Acc_list)
                        Acc_list = []
            if org_flag0 == True:
                break
                
    
    
###################################################################################################################################      
    ##Comment out the specific lookup or general parsing to get an idea of what each one does by itself##    
###################################################################################################################################    
    
    
    #specific lookup of ONLY 1 genbank file
###################################################################################################################################   
################################################################################################################################### 
    #https://docs.python.org/3/tutorial/errors.html
    
    ##Errors/exceptions
#    try:
#        record = db.lookup(gi=556720996) ################This is how you look up a specific record!!!!!!
#        record.format('genbank')
#        print(record.format('genbank'))
#        input('ok?')
#    except IndexError: #If the db.lookup GI value is not in the current database, then you will get an IndexError
#        continue
#    
#    
#    #general parsing
####################################################################################################################################   
####################################################################################################################################        
#    for x in db.values(): #Iterate over DBSeqRecord objects in the namespace [(sub database) aka db_name] 
#            #print(x.annotations)
#            print(x.features)
#            
#            for y in x.features:
#                #print(str(y.qualifiers))
#                print(x.id)
#                print(x.description)
#                #print(x.annotations['])
#                print(str(y.type))
#                print(y.extract(x.seq))
#                print('\n\n')
#                input('ok?')
        
###################################################################################################################################   
################################################################################################################################### 
        