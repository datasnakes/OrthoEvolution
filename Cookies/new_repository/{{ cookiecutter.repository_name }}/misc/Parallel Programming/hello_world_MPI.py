# -*- coding: utf-8 -*-
"""
Created on Thu Sep  8 12:45:57 2016

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
from mpi4py import MPI
from collections import defaultdict
import numpy as np


comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

print(comm)
print(rank)
print(size)
if rank > 0:
    comm.Barrier()
if rank == 0:
    comm.Barrier()
    line_count = 0
    cell_count = 0
    load_list = []

    home = os.getcwd()
    process_size = []
    displacement = []
    load_dict = defaultdict(list)
    
    
    
    loc_load_list = []
    rank_count = 0
    print('check2')
    
    #print('server created: %s' % server)
    os.chdir("/work5/r2294/bin/NCBI_data/vertebrate_mammalian")
    print(os.getcwd())
    
    for file in os.listdir():
        if fnmatch.fnmatch(file, '*rna.gbff'):
            load_list.append(file)
    amt = len(load_list)
    ld_files_rem = amt % (size-2)
    print('remainder: %s' % ld_files_rem)
    ld_files = (amt - ld_files_rem) / (size-2)
    print('ld_files: %s' % ld_files)
    print('load_list: %s' % load_list)
    
    print('amt', amt)
    for x in range(0, size-1):
        print('check')
        print(x)
        print(size-1)
        
        if x == (size-2):
            process_size.append(ld_files_rem)
            displacement.append((ld_files*x))
            print(process_size, displacement)
        else:
            process_size.append(ld_files)
            displacement.append((ld_files*x))
            print(process_size, displacement)
    #tuple(process_size)
    #tuple(displacement)
    
comm.Barrier()  
comm.Bcast(process_size, root=0)

print(process_size)
if rank > 0:  
    my_size = process_size[rank]
    x_loc= np.zeros(my_size, dtype='char')
    comm.Scatterv([load_list, process_size, displacement, MPI.CHAR],x_loc, root = 0)
    print('process_size: %s' % process_size)
    print('displacement: %s' % displacement)
    print('x_local: %s' % x_loc)
    #    comm.Bcast(load_list, root=0)
        
    
    
