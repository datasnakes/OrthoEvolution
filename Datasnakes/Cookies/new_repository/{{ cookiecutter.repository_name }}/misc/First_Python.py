# -*- coding: utf-8 -*-
"""
Created on Wed Aug 17 15:02:37 2016

@author: rgilmore
"""

import sqlite3
from BioSQL import BioSeqDatabase
from Bio import GenBank
import os

server = BioSeqDatabase.open_database(driver = "sqlite3", db = "biosql.db")

db = server.new_database("HTR1A")

dir_list1 = os.listdir()
print(dir_list1)

#for files in  dir_list1:
parser = GenBank.FeatureParser()
iterator = GenBank.Iterator(open("HTR1A_Ailuropoda melanoleuca.gbk"), parser)
db.load(iterator)
db.adaptor.commit()
#input("%s loaded into HTR1A database.  Proceed?")

server.commit()
server.close()


    
