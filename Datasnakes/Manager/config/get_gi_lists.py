"""This script is designed to create a gi list based on the refseq_rna database
for each taxonomy id on the MCSR. It will also convert the gi list into a
binary file which is more efficient to use with NCBI's Standalone Blast tools.
"""
# Modules used
import os
import sys
import platform
import logging as log
from pathlib import Path
from datetime import datetime as d
from time import time
from multiprocessing import Pool
import pandas as pd
from mpi4py import MPI
from Datasnakes.Orthologs.Blast.utils import get_gilists

# Get child process information
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()
machine = platform.node()

# Get args from PBS input
gi_list_path = sys.argv[1]
project_path = sys.argv[2]
index_path = project_path / Path('index')
raw_data_path = project_path / Path('raw_data')

# Set up logging
log.basicConfig(filename=str(raw_data_path / Path("get_gi_lists.log")), level=log.INFO)
log.info("#------------------------------------------------------------------")
log.info("The script name is %s" % os.path.basename(__file__))
log.info("The date and time is currently %s" % str(d.now()))
log.info("#------------------------------------------------------------------")
log.info("Create a GI list for each organism using the taxonomy id and the blastdbcmd tool on the MCSR.")

# Get a taxid list from the file and change to the output directory
tax_ids = pd.read_csv(str(index_path / Path('taxids.csv')), header=None, dtype=str)  # 1st column = tax id's
tax_ids = list(tax_ids[0])

start = time()
with Pool(processes=10) as p:
    p.map(get_gilists, tax_ids)
    log.info("Took {} minutes to get all gi lists.".format((time() - start) / 60))
