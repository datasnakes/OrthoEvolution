# -*- coding: utf-8 -*-
"""
Created on Thu Jan 19 12:58:13 2017

@author: shutchins2

Phylip Test

"""

# List of modules
import os
import pexpect  # I used this to feed input into shell executable


# -----------------------------------------------------------------------------

# Echos all commands in the current shell.
os.system("set -x")

# Create a variable for os.rename
rn = os.rename

# Rename the Phylip input file
rn("")


# -----------------------------------------------------------------------------
# Maximum Likelihood using Phylip executable, dnaml, within unix shell
dnaml = pexpect.spawnu("dnaml infile")
dnaml.sendline("Y\r")
dnaml.waitnoecho()
rn("outfile, _maxlike")
rn("outtree, _maxliketree")

# -----------------------------------------------------------------------------
# Maximum Parsimony using Phylip executable, dnapars, within unix shell
dnapars = pexpect.spawnu("dnapars infile")
dnapars.sendline("Y\r")
dnapars.waitnoecho()
rn("outfile, _maxpars")
rn("outtree, _maxparstree")

# -----------------------------------------------------------------------------
# Distance Matrix using the Phylip executable, dnadist, within unix shell
dnadist = pexpect.spawnu("dnadist infile")
dnadist.sendline("Y\r")
dnadist.waitnoecho()
rn("outfile, _dnadist")