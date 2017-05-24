# -*- coding: utf-8 -*-
"""
Created on Wed Sep 14 10:21:49 2016

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

                 
line_count = 0
org_list = []
gene_list = []
server_db_list = []

server = BioSeqDatabase.open_database(driver="sqlite3", db="/work5/r2294/bin/NCBI_data/vertebrate_mammalian/DB/GPCR_Orthologs_DB.db")
home = os.getcwd()

try:
    #print(server.keys())
    input('is this length right?')
    for db_name in server.keys():
        server_db_list.append(db_name)
    #print(server_db_list)
    print(len(server_db_list))
    input('right?')
    
    with open('/work5/r2294/bin/NCBI_data/Input_Files/Master_Accession_File.csv', 'r') as csvfile:
        Mast_Acc = csv.reader(csvfile, delimiter=',')
        for line in Mast_Acc:
            print('.', end='\r')
            line_count += 1
            cell_count = 0
            if line_count == 1:
                org_count = 0
                for cell in line:
                    org_count += 1
                    if org_count > 1:
                        cell = cell.replace(' ', '_')
                        org_list.append(cell)
                        #print(org_count-2)
                        #print(org_list[org_count-2])
                continue
            for cell in line:
                print('..', end='\r')
                cell_count += 1
                if cell_count == 1:
                    try:
                        os.mkdir('/work5/r2294/bin/NCBI_data/Raw_FASTA_Files/%s' % cell)
                    except FileExistsError:
                        print('file exists')
                    gene_list.append(cell)
                    gene = cell
                    continue
                head, sep, tail = cell.partition('.')
                cell = head
                print(cell)
                Orth_Acc = cell.upper()
                #if cell == '':
#                    if 'Missing_Index.csv' not in os.listdir('/work5/r2294/bin/NCBI_data/Raw_FASTA_Files/'):
#                        with open('/work5/r2294/bin/NCBI_data/Raw_FASTA_Files/Missing_Index.csv', 'a') as csvfile:
#                            Index = csv.writer(csvfile, delimiter=',')
#                            Index.writerow(['GENE', 'ORGANISM', 'GI_NUMBER', 'DB_NAME'])
#                    with open('/work5/r2294/bin/NCBI_data/Raw_FASTA_Files/Missing_Index.csv', 'a') as csvfile:
#                        Index = csv.writer(csvfile, delimiter=',')
#                        Index.writerow([gene, org_list[cell_count-2], 'Missing GI', 'N/A'])
                    #continue
#                for db_name in server.keys():
#                    print('...', end='\r')
                    #print(db_name)
                    
                    #print(db.dbid)
                try:
                    db = server[gene_list[line_count-2]]
                    record = db.lookup(accession=Orth_Acc)
                    #print(Orth_Acc + '\n' + '\n')
                    #print('right gene?')
                    #input(gene_list[line_count-2])
                    #input(org_list[cell_count-2])
                    
                    
                    feature_list = []
                    feat_count = 0
                    for feature in record.features:
                        ref = record.annotations['accessions'][0]
                        GI = record.annotations['gi']
                        feat_count +=1
                        feature_list.append(feature.type)
                        name = str(feature.type)
                        n = str(feature.type)
                        #if feat_count == 1:                         
                            
                            
                       
                        if feature.type in feature_list:
                            x = feature_list.count(feature.type)
                            name = name + str(x)
                        if feature.type == "source":
                            o = feature.qualifiers["organism"]
                            o[0] = '[' + str(o[0]) + ']'
                        if feature.type == "gene":
                            c = 0
                            for x in feature.qualifiers['db_xref']:
                                c += 1
                                if 'GI' in x:
                                    head, sup, tail = x.partition(':')
                        if feature.type == "CDS":
                            
                            with open('/work5/r2294/bin/NCBI_data/Raw_FASTA_Files/%s/%s_%s_%s.ffn' % (gene_list[line_count-2], gene, org_list[cell_count-2], name), 'w') as f:
                                f.write(">" + "gi|" + GI + "|" + "ref" + "|" + ref + '|' + ' ' + record.description + "\n"  + str(feature.extract(record.seq)))
                            with open('/work5/r2294/bin/NCBI_data/Raw_FASTA_Files/%s/%s_%s_%s.faa' % (gene_list[line_count-2], gene, org_list[cell_count-2], 'protein'), 'w') as f:
                                c = 0
                                for x in feature.qualifiers['db_xref']:
                                    c += 1
                                    if 'GI' in x:
                                        head, sup, tail = x.partition(':')
                                f.write(">" + "gi|" + tail + "|" + "ref" + "|" +  str(feature.qualifiers['protein_id'][0]) + '| ' + str(feature.qualifiers['product'][0]) + ' ' + o[0] + '\n' + str(feature.qualifiers['translation'][0]))
                            
                            
                            with open('/work5/r2294/bin/NCBI_data/Raw_FASTA_Files/%s/MASTER_%s_%s.ffn' % (gene_list[line_count-2], gene, n), 'a') as f:
                                f.write(">" + "gi|" + GI + "|" + "ref" + "|" + ref + '|' + ' ' + record.description + "\n"  + str(feature.extract(record.seq))+ "\n\n")
                            
                            with open('/work5/r2294/bin/NCBI_data/Raw_FASTA_Files/%s/MASTER_%s_%s.faa' % (gene_list[line_count-2], gene, 'protein'), 'a') as f:
                                c = 0
                                for x in feature.qualifiers['db_xref']:
                                    c += 1
                                    if 'GI' in x:
                                        head, sup, tail = x.partition(':')
                                f.write(">" + "gi|" + tail + "|" + "ref" + "|" +  str(feature.qualifiers['protein_id'][0]) + '| ' + str(feature.qualifiers['product'][0]) + ' ' + o[0] + '\n' + str(feature.qualifiers['translation'][0]) + '\n\n')

                            print(SeqIO.read('/work5/r2294/bin/NCBI_data/Raw_FASTA_Files/%s/%s_%s_%s.ffn' % (gene_list[line_count-2], gene, org_list[cell_count-2], name), 'fasta'))
                            
                        elif feature.type == "misc_feature":
                            
                            with open('/work5/r2294/bin/NCBI_data/Raw_FASTA_Files/%s/%s_%s_%s.ffa' % (gene_list[line_count-2], gene, org_list[cell_count-2], name), 'w') as f:
                                f.write(">" + "gi|" + GI + "|" + "ref" + "|" + ref + '|' + ' ' + record.description + "\n" + 'Feature: ' + str(feature.qualifiers['note'][0]) + '\n'  + str(feature.extract(record.seq)) + "\n")
                            
                            with open('/work5/r2294/bin/NCBI_data/Raw_FASTA_Files/%s/MASTER_%s_%s.ffa' % (gene_list[line_count-2], gene, n), 'a') as f:
                                f.write(">" + "gi|" + GI + "|" + "ref" + "|" + ref + '|' + ' ' + record.description + "\n"  + 'Feature: ' + str(feature.qualifiers['note'][0]) + '\n' + str(feature.extract(record.seq)) + "\n\n")
                            print(SeqIO.read('/work5/r2294/bin/NCBI_data/Raw_FASTA_Files/%s/%s_%s_%s.ffa' % (gene_list[line_count-2], gene, org_list[cell_count-2], name), 'fasta'))
                            
                        elif feature.type != "variation":
                            with open('/work5/r2294/bin/NCBI_data/Raw_FASTA_Files/%s/%s_%s_%s.fasta' % (gene_list[line_count-2], gene, org_list[cell_count-2], name), 'w') as f:
                                f.write(">" + "gi|" + GI + "|" + "ref" + "|" + ref + '|' + ' ' + record.description + "\n"  + str(feature.extract(record.seq)) + "\n")
                            
                            with open('/work5/r2294/bin/NCBI_data/Raw_FASTA_Files/%s/MASTER_%s_%s.fasta' % (gene_list[line_count-2], gene, n), 'a') as f:
                                f.write(">" + "gi|" + GI + "|" + "ref" + "|" + ref + '|' + ' ' + record.description + "\n"  + str(feature.extract(record.seq)) + "\n\n")
                            print(SeqIO.read('/work5/r2294/bin/NCBI_data/Raw_FASTA_Files/%s/%s_%s_%s.fasta' % (gene_list[line_count-2], gene, org_list[cell_count-2], name), 'fasta'))
                        
                        
                            
#                    with open(('/work5/r2294/bin/NCBI_data/Raw_FASTA_Files/%s/%s_%s.fasta' % (gene_list[line_count-2], gene, org_list[cell_count-2])), 'w') as GB_File:
#                        GB_File.write(record.format('fasta'))
#                        print(record.format('fasta'), 'saved')
                    
                    print(os.listdir('/work5/r2294/bin/NCBI_data/Raw_FASTA_Files/%s' % gene_list[line_count-2]))
                    #print('right gene?')
                    #input(gene_list[line_count-2])
#                        if 'Index.csv' not in os.listdir('/work5/r2294/bin/NCBI_data/Raw_FASTA_Files/'):
#                            with open('/work5/r2294/bin/NCBI_data/Raw_FASTA_Files/Index.csv', 'a') as csvfile:
#                                Index = csv.writer(csvfile, delimiter=',')
#                                Index.writerow(['GENE', 'ORGANISM', 'GI_NUMBER', 'DB_NAME'])
#                        with open('/work5/r2294/bin/NCBI_data/Raw_FASTA_Files/Index.csv', 'a') as csvfile:
#                            Index = csv.writer(csvfile, delimiter=',')
#                            Index.writerow([gene, org_list[cell_count-2], Orth_Acc, db_name])
                                    
                                    
                    #input('is this ok?')
                    print('....', end='\r')
                    #break
                except IndexError:
                    #print('cell', cell, 'gene list', gene_list[line_count - 2], 'org list', org_list[cell_count - 2], 'cell count', cell_count)
                    #print('index error')
                    #input('This gene eneds to be updated')
                    print('index error')
except:
    server.rollback()
    raise