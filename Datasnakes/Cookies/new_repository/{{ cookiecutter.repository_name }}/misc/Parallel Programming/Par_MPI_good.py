# -*- coding: utf-8 -*-
"""
Created on Tue Sep 13 16:36:20 2016

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
if rank == 0: ############################    RANK 0
    sub_list = []
    load_list = []
    os.chdir("/work5/r2294/bin/NCBI_data/vertebrate_mammalian")
    print(os.getcwd())
    
    for file in os.listdir():
        if fnmatch.fnmatch(file, '*rna.gbff'):
            load_list.append(file)
    #print(load_list)
    sub_list = load_list
   # print('sub_list', sub_list)
     
    for x in range(1, size):
        print('x: ', x)
        print(sub_list[0], 'root')
        comm.send(sub_list[0], dest=x, tag = 13)
        sub_list.pop(0)
    
elif size > rank > 0: ############################    RANK 1
    db = comm.recv(source = 0, tag = 13)
    print('db: ', db, 'rank: ', rank)
    sub_list = None
sub_list = comm.bcast(sub_list, root=0)
    #Do work Here

    
#    send_flag = True
#    print('check', rank)
#    comm.send(send_flag, dest=0, tag = 44)
#    print('check', rank)
#    comm.send(rank, dest=0, tag = 33)
#    db = comm.recv(source = 0, tag = 55)

  
while sub_list != []:
    if rank == 0:
        req = comm.irecv(source = MPI.ANY_SOURCE, tag = 44)
        print('check', rank)
        send_flag = req.wait()
        if send_flag == True:
            source = comm.recv(source= MPI.ANY_SOURCE, tag = 33)
            comm.send(sub_list[0], dest=source, tag = 55)
            print(sub_list[0])
            sub_list.pop(0)
            send_flag = False
    if size > rank > 0:
        send_flag = True
        print('check', rank)
        comm.send(send_flag, dest=0, tag =44)
        print('check', rank)
        comm.send(rank, dest=0, tag=33)
            

#bcast(sub_list, root = 0)
            