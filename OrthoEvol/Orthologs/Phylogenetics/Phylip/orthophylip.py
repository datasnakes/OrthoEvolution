import os
import pexpect  # I used this to feed input into shell executable
import sys


class Phylip(object):
    def __init__(self, inputfile):
    """The input file should be a phylip formatted multiple sequence
    alignment."""

        self._rename = os.rename
        if sys.platform == 'win32' or 'win64':
            sys.exit("This module is strictly for use on Linux at the moment.")

        self.inputfile = inputfile

        # Rename the input file to infile
        self._rename(self.inputfile, "infile")
        self.inputfile = "infile"

    def dnapars(self, outfile, outtree):
        """Maximum Parsimony using Phylip executable, dnapars,
within unix shell.

        :param outfile: 
        :param outtree: 

        """
        dnapars = pexpect.spawnu("dnapars infile")
        dnapars.sendline("Y\r")
        dnapars.waitnoecho()
        self._rename("outfile", outfile + "_dnapars_output")
        self._rename("outtree", outtree + "_maxparsimony_tree")

    def dnaml(self, outfile, outtree):
        """Maximum Likelihood using dnaml within a unix shell.

        :param outfile: 
        :param outtree: 

        """
        dnaml = pexpect.spawnu("dnaml infile")
        dnaml.sendline("Y\r")
        dnaml.waitnoecho()
        self._rename("outfile", outfile + "_dnaml_output")
        self._rename("outtree", outtree + "_maxlikelihood_tree")

    def dnadist(self, dnadist_output):
        """

        :param dnadist_output: 

        """
        dnadist = pexpect.spawnu("dnadist infile")
        dnadist.sendline("Y\r")
        dnadist.waitnoecho()
        self._rename("outfile", dnadist_output + "_dnadist")
