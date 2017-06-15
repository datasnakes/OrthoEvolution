# -*- coding: utf-8 -*-
"""
Created on Fri Sep  2 17:35:03 2016

@author: rgilmore
"""

from subprocess import Popen
import subprocess
import time

# If you want to be emailed by the system, include these in job_string:
#PBS -M your_email@address
#PBS -m abe # (a = abort, b = begin, e = end)

for i in range (1, 3):
    
    output, inp = subprocess.run('qsub', stdout=subprocess.PIPE)
    job_name = "Refseq_Full_Load_%d" % i
    walltime = "10:00:00"
    processors = "nodes=2:ppn=2"
    command = "./DB_creation.py -n %d" % i
    print('Subprocess%d has started' % i)
    
    job_string = """#!/bin/bash
    #PBS -M rgilmore@umc.edu
    #PBS -m abe # (a = abort, b = begin, e = end)
    #PBS -N %s
    #PBS -l walltime = %s
    #PBS -l %s
    #PBS -o ./output/%s.err
    #PBS -e ./error/%s.err
    cd $PBS_O_WORKDIR
    %s""" % (job_name, walltime, processors, job_name, job_name, command)
    
    inp.write(job_string)
    inp.close()
    
    print(job_string)
    print(output.read())
    
    time.sleep(0.1)
    