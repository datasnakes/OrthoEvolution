# -*- coding: utf-8 -*-
"""
Created on Mon Sep 12 12:14:06 2016

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
import numpy


comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

print(comm)
print(rank)
print(size)
if rank > 0:
    comm.Barrier()
if rank == 0:
    
    line_count = 0
    cell_count = 0
    load_list = []
    #server = BioSeqDatabase.open_database(driver="sqlite3", db="TEST_DB.db")
    home = os.getcwd()
    process_size = []
    displacement = []
    committed_db = []
    loaded_list = []
    current_load = ''
    load_dict = defaultdict(list)
    
    
    
    loc_load_list = []
    rank_count = 0

#def load_check(db_name, rank, load_dict):
#    for x in range(1, size):
#        if db_name in load_dict[x]:
#            return True
#    return False
            
if rank == 0:
    #print('server created: %s' % server)
    os.chdir("/work5/r2294/bin/NCBI_data/vertebrate_mammalian")
    print(os.getcwd())
    
    for file in os.listdir():
        print(file)
        if fnmatch.fnmatch(file, '*rna.gbff'):
            load_list.append(file)
    amt = len(load_list)
    ld_files_rem = amt % (size-2)
    print('remainder: %s' % ld_files_rem)
    ld_files = (amt - ld_files_rem) / (size-2)
    print('ld_files: %s' % ld_files)
    #print('load_list: %s' % load_list)
    
   # print('amt', amt)
    for x in range(1,size+1):
        print('check')
        print(x)
        print(size-1)
        if x == (size-1):
            process_size.append(ld_files_rem)
            displacement.append((ld_files*x)-1)
            print(process_size, displacement)
        else:
            process_size.append(ld_files)
            displacement.append((ld_files*x)-1)
            print(process_size, displacement)
    x_local = numpy.zeros(size-1)
    tuple(process_size)
    tuple(displacement)
    comm.Scatterv([load_list, process_size, displacement, MPI.DOUBLE], x_local, root = 1)
    print('process_size: %s' % process_size)
    print('displacement: %s' % displacement)
    print('x_local: %s' % x_local)
    comm.Bcast(load_list, root=0)
        
    
    comm.Barrier()

#if rank > 0:
    #comm.Barrier()

    #try:
#        loc_load_list = load_list
#        while rank > rank_count:
#            rank_count += 1
#            loc_load_list.pop(0)
        
        
            
#        for file in loc_load_list:
#            
#            check = load_check(file, rank, load_dict)
#            db_name = file
#            db_flag = False
#            
#            if file in server.keys():
#                continue
#            else:
#                if check is False:
#                    load_dict[rank].append(db_name)
#                    comm.Bcast(load_dict)
#                try:
#                    
#            current_load = db_name            
#            
##            if current_load in loaded_list:
##                continue
##            current_load = db_name
##            comm.Bcast(current_load, root=0)
##            loaded_list = loaded_list.append(current_load)
    #except:
     #       server.rollback()
      #      raise