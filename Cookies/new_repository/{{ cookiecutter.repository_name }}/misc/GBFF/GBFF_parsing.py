# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.

"""

from Bio import SeqIO
import fnmatch
from BioSQL import BioSeqDatabase


#server = BioSeqDatabase.open_database(driver="sqlite3", db ="Test_DB.db")
#
#
#
#print('server created: %s' % server)
#
##server.remove_database("GB_files")
#del server["GB_files"]
#server.commit()
#db = server.new_database("GB_files")
##loading = BioSeqDatabase.Loader.DatabaseLoader(server, db)
#
#
#print('database created GB_files: %s' % db)
#server.commit()
#for x in server.keys():
#    print(x)



output1 = open("SeqIO_Class.txt", "w")
output2 = open("SeqIO_modules.txt", "w")
temp = []
file = "vertebrate_mammalian.535.rna.gbff"
rec = 0
count = 0
#for seq_record in SeqIO.parse(file, "genbank"):
#print(seq_record.id.count("."))
#coun = db.load(SeqIO.parse(file, "genbank"))


#input('this is what the genbank file looks liek')
#print(coun)
#print("this database contains %i records" % len(db))
#input('this is what the database file looks like')
#server.commit()
    
input('next')
for record in SeqIO.parse(file, "genbank"):
    print(dir(record))
    output1.write("Here is what's available in the record class\n\n")
    
    for x in dir(record):
        if fnmatch.fnmatch(x, '_*_') == True:
            continue
        output1.write('\t' + str(x) + '\n')
        output1.write(("\t\tHere's what's available in the \'%s\' module" % x) + '\n')
        count = 0
    output2.write("Here is what's available in the record class\n\n")
    
    for x in dir(record):
        for y in dir(x):
            count = count + 1
            if fnmatch.fnmatch(y, '_*_') == True:
                continue
            #if count == 5:
                #print('')
            if count > 1:
                output2.write('\t\t\t')
            output2.write(str(y))
            
    print(record.annotations.keys())
    print(record.annotations.values())
    print(record.seq)
    print(dir(record.seq))
    rec = rec + 1
    print("this is record %s" % str(rec))
    print(record)
    input('')
    print(record.annotations['gi'])
    ct = 0
    for this in record.features:
        print(this.type)
        print(this.qualifiers)
        #print(this.type.qualifiers)
        print(this)
        input("ok")
    for this in record.features:
        if this.type == 'gene':
            ct = ct + 1
            name = str(this.qualifiers['gene'])
            print(name)
            
            print("\n This is feature #: " + str(ct))
            print("And this is the type: " + str(this.type))
            #print(this.qualifiers['note'])
            print(this.extract(record.seq))
            print("\n")
