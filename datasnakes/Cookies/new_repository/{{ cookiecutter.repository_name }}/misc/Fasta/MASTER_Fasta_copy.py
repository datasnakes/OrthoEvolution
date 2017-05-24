# -*- coding: utf-8 -*-
"""
Created on Tue Sep 20 10:43:44 2016

@author: rgilmore
"""

import os
import fnmatch
import shutil
import csv

home = os.getcwd()
folder_list = []

os.chdir('/work5/r2294/bin/NCBI_data/')
if os.path.isdir('/work5/r2294/bin/NCBI_data/Raw_FASTA_Files') == False:
    os.mkdir('/work5/r2294/bin/NCBI_data/Raw_FASTA_Files')
os.chdir('/work5/r2294/bin/NCBI_data/Raw_FASTA_Files')
with open('/work5/r2294/Database/Gene_List.csv', 'r')as csvfile:
    gene = csv.reader(csvfile, delimiter=',')
    for line in gene:
        for cell in line:
            folder_list.append(cell)
    print(folder_list)
    input('Is this right?')

#Set up Directories for MASTER FASTA FILES
home = '/work5/r2294/bin/NCBI_data/Raw_FASTA_Files/'
Mast_home = '/work5/r2294/bin/NCBI_data/Raw_FASTA_Files/MASTER_FASTA_FILES/'
if os.path.isdir(Mast_home) == False:
    os.mkdir(Mast_home)
if os.path.isdir(Mast_home + 'MASTER_CDS') == False:
    os.mkdir(Mast_home + 'MASTER_CDS')
if os.path.isdir(Mast_home + 'MASTER_Protein') == False:
    os.mkdir(Mast_home + 'MASTER_Protein')
if os.path.isdir(Mast_home + 'MASTER_Misc_Features') == False:
    os.mkdir(Mast_home + 'MASTER_Misc_Features')
if os.path.isdir(Mast_home + 'MASTER_Other') == False:
    os.mkdir(Mast_home + 'MASTER_Other')
if os.path.isdir(Mast_home + 'MASTER') == False:
    os.mkdir(Mast_home + 'MASTER')


for folder1 in os.listdir(Mast_home):
    if folder1 == 'MASTER':
        continue
    if folder1 == 'MASTER_Other':
        continue
    os.chdir(Mast_home)
    x = str(folder1)
    x = x.replace('MASTER_', '')
    folder2 = x
    print(folder2)
    if x == 'Misc_Features':
        x = x[:-1]
        folder2 = x
        print(folder2)
    if os.path.isdir(Mast_home + folder1 + '/' + folder2) == False:
        os.mkdir(Mast_home + folder1 + '/' + folder2)
#Creates directories first if they haven't already been created
    for folder3 in folder_list:
        if os.path.isdir(Mast_home + folder1 + '/' + folder2 + '/' + folder3) == False:
            os.mkdir(Mast_home + folder1 + '/' + folder2 + '/' + folder3)
    
#Finds less important sections (gene, source, STS, exon, etc.) and puts them in an Other folder
    if 'Other' in folder1:
        for folder3 in folder_list:            
            for file in os.listdir(home + folder3):
                os.chdir(home + folder3)
                
                if fnmatch.fnmatch(file, '*misc_feature*') == False:
                    if fnmatch.fnmatch(file, '*CDS*') == False:
                        if fnmatch.fnmatch(file, '*protein*') == False:
                            shutil.copy2(file, (Mast_home + folder1 + '/' + folder2 + '/' + folder3))
                            print(folder1, '%s files moved' % folder3)
                            if fnmatch.fnmatch(file.lower(), ('master*')) == True:
                                shutil.copy2(file, (Mast_home + folder1))
        continue
    
    for folder3 in folder_list:        
        for file in os.listdir(home + folder3):
            os.chdir(home + folder3)
#Matches CDS, Protein, and Misc_Features FASTA FILES and puts them into corresponding folders
            if fnmatch.fnmatch(file.lower(), ('*%s*' % folder2.lower())) == True:
                shutil.copy2(file, (Mast_home + folder1 + '/' + folder2 + '/' + folder3))
                print(folder1, ' %s files moved' % folder3)
                if fnmatch.fnmatch(file.lower(), ('master*')) == True:
                    shutil.copy2(file, (Mast_home + folder1))