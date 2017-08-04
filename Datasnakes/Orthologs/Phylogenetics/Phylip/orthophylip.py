import os
import pexpect  # I used this to feed input into shell executable
import sys

# Create a variable for os.rename
rn = os.rename


class Phylip(object):
    def __init__(inputfile):
        """ The input file should be a phylip formatted multiple sequence
        alignment."""
        if sys.platform == 'win32' or 'win64':
            sys.exit("This module is strictly for use on Linux at the moment.")

        # Rename the input file to infile
        rn(inputfile, "infile")

    def dnapars(gene):
        """ Maximum Parsimony using Phylip executable, dnapars,
        within unix shell."""
        dnapars = pexpect.spawnu("dnapars infile")
        dnapars.sendline("Y\r")
        dnapars.waitnoecho()
        rn("outfile", gene + "_maxpars")
        rn("outtree", gene + "_maxparstree")

    def dnaml(gene):
        """Maximum Likelihood using Phylip executable, dnaml, within a unix shell. """
        dnaml = pexpect.spawnu("dnaml infile")
        dnaml.sendline("Y\r")
        dnaml.waitnoecho()
        rn("outfile", gene + "_maxlike")
        rn("outtree", gene + "_maxliketree")

    def dnadist(gene):
        dnadist = pexpect.spawnu("dnadist infile")
        dnadist.sendline("Y\r")
        dnadist.waitnoecho()
        rn("outfile", gene + "_dnadist")
