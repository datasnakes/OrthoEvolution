# -*- coding: utf-8 -*-
"""
Created on Thu Sep  1 11:28:04 2016

@author: rgilmore
"""

import csv
import os

MAF = "Master_Accession_File.csv"
MGF = "Master_GI_File.csv"

line_count = 0
cell_count = 0
Gene = ''
Species = ''
Accession = ''
Org = {Species:Accession}
Full_dict = {Gene:Org}

with open(MAF, 'r') as csvAC, open(MGF, 'r') as csvGI:
    Acc_File = csv.reader(csvAC)
    GI_File = csv.reader(csvGI)
    
    for (A_line, G_line) in (Acc_File, GI_File):
        line_count += 1
        for(A_cell, G_cell) in (A_line, G_line):
            cell_count += 1
            
            if line_count == 1 and cell_count == 1:
                continue
            elif line_count == 1 and cell_count > 1:
                org_list.append(A_cell)
            elif line_count > 1 and cell_count == 1:
                gene_list.append(A_cell
            elif line_count > 1 and cell_count > 1:
                
            cell_count = 0
