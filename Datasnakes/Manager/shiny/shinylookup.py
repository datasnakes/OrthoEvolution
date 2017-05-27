
"""

"""
# Libraries:

import os
import sys
from pathlib import Path

from Bio import SeqIO
from BioSQL import BioSeqDatabase

from dir_mana import dir_mana
from lister import Lister

# Custom Class Initializations
# :
# Use directory_management() class here so that we can stay organized
# and more easily access the proper directories on command
home = os.getcwd()
project = "GPCR-Orthologs-Project"
user = "Work"
where = dir_mana(home, project)
# Use lister() class here so that we can easily access our Master RNA
# Accession File
what = Lister('MAFV3.1.csv')  # Always make sure this file name is correct

# Add a path that contains custom libraries for import
# os.sys.path.append()
# Gene and Organism choices from the R script are used here for database lookup
Choices = sys.argv

Accession = what.gene_dict[Choices[1]][Choices[2]]
Accession, Sup, Version = (str(Accession)).partition('.')
Accession = Accession.upper()

GenBank = '%s_%s.gbk' % (Choices[1], Choices[2])
Attributes = '%s_%s.kv' % (Choices[1], Choices[2])
Features = '%s_%s.feat' % (Choices[1], Choices[2])
Annotations = '%s_%s.anno' % (Choices[1], Choices[2])
Comment = '%s_%s.comm' % (Choices[1], Choices[2])
Reference = '%s_%s.ref' % (Choices[1], Choices[2])

db_list = []
for file in os.listdir(where.APP_DATA):
    if '.db' in str(file):
        db_list.append(file)
for db_name in db_list:
    server = BioSeqDatabase.open_database(
        driver='sqlite3', db=str(
            where.APP_DATA / Path(db_name)))
    for sub_db_name in server.keys():
        db = server[sub_db_name]
        try:
            record = db.lookup(accession=Accession)
            with open(where.APP_DATA / Path(GenBank), 'w') as GB_file:
                GB_file.write(record.format('genbank'))
                print(GB_file.name, 'created')
            server_flag = True
        except IndexError:
            continue

record = SeqIO.parse(where.APP_DATA / Path(GenBank), 'genbank')
Choice_list = ''
Display_list = ''
item_list = []
type_count = 0
if os.path.isfile(where.APP_DATA / Path(Attributes)):
    os.remove(where.APP_DATA / Path(Attributes))

for record in SeqIO.parse(str(where.APP_DATA / Path(GenBank)), 'genbank'):
    with open(where.APP_DATA / Path(Attributes), 'w', newline='\n') as file:
        for key, value in record.__dict__.items():
            key = str(key)
            value = str(value)
            file.write('%s#%s\n' % (key, value))
    with open(where.APP_DATA / Path(Features), 'w', newline='\n') as file:
        for item in record.features:
            if item.type == 'variation':
                continue
            if item.type in item_list:
                type_count += 1
            item_list.append(item.type)
            file.write(
                '%s#%s\n' %
                (item.type +
                 str(type_count),
                    item.extract(
                     record.seq)))
    with open(where.APP_DATA / Path(Annotations), 'w', newline='\n') as file:
        for k, v in record.annotations.items():
            if k == 'comment' or k == 'references':
                continue
            file.write('%s#%s\n' % (k, v))
            # print(k, v)
try:
    with open(where.APP_DATA / Path(Reference), 'w') as file:
        for item in record.annotations['reference']:
            file.write(item.pubmed_id)
except KeyError:
    with open(where.APP_DATA / Path(Reference), 'w') as file:
        file.write('No references in the GenBankFile')
with open(where.APP_DATA / Path(Comment), 'w') as file:
    try:
        file.write(record.annotations['comment'])
    except KeyError:
        file.write('No comment in the GenBank File')
