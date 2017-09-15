# -*- coding: utf-8 -*-
"""
Date created: Tue Mar  7 15:38:34 2017
Author: S. Hutchins

Script description: Testing different Clustal Omega parameters.

"""
import os
from ClustalParams import COP1, COP2, COP3, COP4, COP5, COP6, COP7

# Assess current working directories and establishing home and output directories.

a = '/work5/r2295/bin/ClustalParamsTest/' # Home directory
home = a
os.chdir(home) # Directory change to output directory

paramnames = ["params1", "params2", "params3", "params4", "params5", "params6", "params7"]
parameters = [COP1, COP2, COP3, COP4, COP5, COP6, COP7]

for name, param in zip(paramnames, parameters):
    # Create directories for alignment files.
    b = a + name + "-output"
    os.makedirs('%s' % b, exist_ok=True) # Create a directory or don't if it exists.

    # Go to the CDS directory
    c = a + "InputCDS"

    os.chdir(c) # Change to cds.fasta file directory
    os.listdir(c) # Make a list of the files in the current directory
    print("➜ Current CDS/Gene directory: "+ os.getcwd() + "\n")  # Print current working directory

    # Echos all commands in the current shell.
    os.system("set -x")

    # Copies the profile.fasta and concatenated cds.fasta file to output dir.
    os.system("cp {APOL1_cds_nucl.fasta,profile.fasta} " + b + "/")

    # Directory change to output directory
    os.chdir(b)

    # Run Clustal Omega Commandline tool
    param(in_file="APOL1_cds_nucl.fasta", out_file=name + ".fasta", outfmt="fasta", logfile=name + ".log")