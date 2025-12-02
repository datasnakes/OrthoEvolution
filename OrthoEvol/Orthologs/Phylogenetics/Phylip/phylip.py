import os
import sys
import shutil

import pexpect  # I used this to feed input into shell executable
from Bio import AlignIO

from OrthoEvol.Tools.logit import LogIt


class Phylip(object):
    """A class that serves as a wrapper for the Phylip executable."""

    def __init__(self, infile):
        """Initialize the Phylip class.

        :param infile: A phylip formatted multiple sequence alignment.
        :type infile: str
        :raises OSError: If the operating system is not Linux.
        :raises ValueError: If the input file is not in valid phylip format.
        """
        # Set up logging first
        self.phylip_log = LogIt().default(logname="Phylip", logfile=None)
        # Raise error if OS is not linux
        if sys.platform != 'linux':
            err_msg = "This module is strictly for use on Linux at the moment."
            raise OSError(err_msg)
        # Validate format and set infile
        if not self._validate_format(infile):
            raise ValueError(f"Invalid phylip format for file: {infile}")
        self.infile = infile
        self._rename = os.rename

    def _validate_format(self, infile):
        """Validate the format of the Phylip file

        :param infile: A phylip formatted multiple sequence alignment.
        :type infile: str
        """
        try:
            with open(infile, 'r') as f:
                AlignIO.read(f, "phylip")
            return True
        except (ValueError, IOError) as e:
            self.phylip_log.exception(e)
            return False

    def _temp_infile(self, infile):
        """Create a temporary infile named infile.

        :param infile:  A phylip formatted multiple sequence alignment.
        """
        shutil.copyfile(infile, "infile")
        temp_infile = "infile"
        return temp_infile

    def dnapars(self, outfile, outtree):
        """Generate a maximum parsimony tree using dnapars.

        :param outfile: Standard output filename.
        :type outfile: str
        :param outtree: Name of maximum parsimony tree.
        :type outtree: str
        """
        infile = self._temp_infile(infile=self.infile)
        try:
            dnapars = pexpect.spawnu("dnapars %s" % infile)
            dnapars.sendline("Y\r")
            dnapars.waitnoecho()
        except pexpect.EOF as e:
            self.phylip_log.error(dnapars.read())
            self.phylip_log.exception(e)
        else:
            self.phylip_log.info(dnapars.read())
            self._rename("outfile", outfile)
            self._rename("outtree", outtree)
        finally:
            os.remove(infile)

    def dnaml(self, outfile, outtree):
        """Generate a maximum likelihood tree using dnaml.

        :param outfile: Standard output filename.
        :type outfile: str
        :param outtree: Name of maximum likelihood tree.
        :type outtree: str
        """
        infile = self._temp_infile(infile=self.infile)
        try:
            dnaml = pexpect.spawnu("dnaml %s" % infile)
            dnaml.sendline("Y\r")
            dnaml.waitnoecho()
        except pexpect.EOF as e:
            self.phylip_log.error(dnaml.read())
            self.phylip_log.exception(e)
        else:
            self.phylip_log.info(dnaml.read())
            self._rename("outfile", outfile)
            self._rename("outtree", outtree)
        finally:
            os.remove(infile)

    def dnadist(self, outfile):
        """Generate a distance matrix using dnadist.

        :param outfile: Distance matrix output filename.
        :type outfile: str
        """
        infile = self._temp_infile(infile=self.infile)
        try:
            dnadist = pexpect.spawnu("dnadist %s" % infile)
            dnadist.sendline("Y\r")
            dnadist.waitnoecho()
        except pexpect.EOF as e:
            self.phylip_log.error(dnadist.read())
            self.phylip_log.exception(e)
        else:
            self.phylip_log.info(dnadist.read())
            self._rename("outfile", outfile)
        finally:
            os.remove(infile)
