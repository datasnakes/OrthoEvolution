##############################################################################
# PyCharm Community Edition
# -*- coding: utf-8 -*-
"""
GPCR_Orthologs
db2gbk.py updated on 1/3/2017 at 1:03 PM
##############################################################################

    Input:  Master_Accession_File in a specified format (described in a seperate
    text file) and the proper database files.

    Output:  Extract target GenBank files from the database files that were
    created with ftp2db

    Description:

##############################################################################
@author: rgilmore
"""
##############################################################################
# Libraries:

import os
import shutil

from Bio import SeqIO
from BioSQL import BioSeqDatabase

from dir_mana import dir_mana
from lister import Lister

# //TODO-ROB Add a progress bar to the pipeline
##############################################################################
# Custom Class Initializations
# :
# Use directory_management() class here so that we can stay organized
# and more easily access the proper directories on command
home = os.getcwd()
project = "GPCR-Orthologs-Project"
user = "rgilmore"
where = dir_mana(home, project)
# Use lister() class here so that we can easily access our Master RNA
# Accession File


# Add a path that contains custom libraries for import
# os.sys.path.append()
##############################################################################
# Global Initializations:

##############################################################################
# //TODO-ROB make code more versatile for multiple projects or even single queries
class Db2Gbk(object):

    def __init__(self, m_file='', path='/', genbank_update=False):
        self.M_file = m_file
        # Always make sure this file name is correct
        self.what = Lister('MAFV3.1.csv')
        self.GenBank_update = genbank_update

        path_list = where.path_list_make(where.VERT_MAM, where.NCBI_DATA)
        print('path_list: ', path_list)
        if path == '/':
            self.path, t = self.db2gbk_mkdir(
                where.VALLENDER_DATA, path_list, self.GenBank_update)
        else:
            self.path = path
            print('self.path: ', self.path)

    def gbk_gather(self):
        os.chdir(where.VERT_MAM + '/Databases')
        db_name = []
        for file in os.listdir(os.getcwd()):
            if file.endswith('.db'):
                db_name.append(str(file))

        for G_key, G_value in self.what.tier_frame_dict.items():
            Tier = G_key
            os.chdir(self.path)
            os.mkdir(Tier)
            os.chdir(Tier)
            Tier_path = os.getcwd()
            for Gene in self.what.tier_frame_dict[Tier].T:
                os.chdir(Tier_path)
                os.mkdir(Gene)
                os.chdir(Gene)
                for Organism in self.what.org_list:
                    Accession = str(self.what.gene_dict[Gene][Organism])
                    Accession, Sup, Version = Accession.partition('.')
                    Accession = Accession.upper()
                    server_flag = False
                    for name in db_name:
                        if server_flag is True:
                            break
                        name = str(name)
                        server = BioSeqDatabase.open_database(
                            driver='sqlite3', db=where.VERT_MAM + ('/Databases/%s' % name))
                        for sub_db_name in server.keys():
                            db = server[sub_db_name]

                            try:
                                record = db.lookup(accession=Accession)
                                with open('%s_%s.gbk' % (Gene, Organism), 'w') as GB_file:
                                    GB_file.write(record.format('genbank'))
                                    print(GB_file.name, 'created')
                                server_flag = True
                                break
                            except IndexError:
                                print('Index Error')
                                continue

    def gbk_upload(self):
        t_count = 0
        os.chdir(self.path)
        print(os.getcwd())
        if os.path.isdir(self.path + '/Databases') is False:
            os.mkdir('Databases')
        for tier in os.listdir(os.getcwd()):
            if tier == 'Databases':
                continue
            db_name = str(tier) + '.db'
            if os.path.isfile(self.path + '/Databases/' + db_name) is False:
                print('Copying Template BioSQL Database...  '
                      'This may take a few minutes...')
                shutil.copy2(where.Templates + '/Template_BioSQL_DB.db',
                             self.path + '/Databases/%s' % db_name)
            else:
                os.remove(self.path + '/Databases/' + db_name)
                print('Copying Template BioSQL Database...  '
                      'This may take a few minutes...')
                shutil.copy2(where.Templates + '/Template_BioSQL_DB.db',
                             self.path + '/Databases/%s' % db_name)

            server = BioSeqDatabase.open_database(
                driver='sqlite3', db=(
                    self.path + '/Databases/' + db_name))
            os.chdir(tier)
            for gene in os.listdir(os.getcwd()):
                os.chdir(gene)
                sub_db_name = gene
                for file in os.listdir(os.getcwd()):
                    try:
                        if sub_db_name not in server.keys():
                            server.new_database(sub_db_name)
                        db = server[sub_db_name]
                        count = db.load(SeqIO.parse(file, 'genbank'))
                        server.commit()
                        print('Server Commited %s' % sub_db_name)
                        print('%s database loaded with %s.' % (db.dbid, file))
                        print(
                            "That file contains %s genbank records." %
                            str(count))
                        t_count = t_count + count
                        print(
                            'The total number of files loaded so far is %i.' %
                            t_count)
                    except BaseException:
                        server.rollback()
                        try:
                            del server[sub_db_name]
                            server.commit()
                        except BaseException:
                            raise
                        raise
                os.chdir('..')
            os.chdir('..')

    def db2gbk_mkdir(self, path, p_list, update):
        if update is True:
            path = where.dir_archive(path, p_list)
        else:
            path = where.dir_make(path, p_list)
        return path
