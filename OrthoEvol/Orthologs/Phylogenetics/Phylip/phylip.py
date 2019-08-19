import os
import sys
import shutil

import pexpect  # I used this to feed input into shell executable

from OrthoEvol.Tools.logit import LogIt


class Phylip(object):
    """A class that serves as a wrapper for the Phylip excecutable."""

    def __init__(self, infile):
        """Initialize the Phylip class.

        :param infile: A phylip formatted multiple sequence alignment.
        """
        self.infile = infile
        self._rename = os.rename
        # Set up logging
        self.phylip_log = LogIt().default(logname="Phylip", logfile=None)
        # Raise error is OS is not linux
        if sys.platform != 'linux':
            err_msg = "This module is strictly for use on Linux at the moment."
            raise OSError(err_msg)

    def _validate_format(self, infile):
        """Validate the format of the Phylip file

        :param infile: A phylip formatted multiple sequence alignment.
        :type infile: str
        """
        pass

    def _temp_infile(self, infile):
        """Create a temporary infile named infile.

        :param infile:  A phylip formatted multiple sequence alignment.
        """
        shutil.copyfile(self.infile, "infile")
        temp_infile = "infile"
        return temp_infile

    def dnapars(self, outfile, outtree):
        """Generate a maximum parsimony tree using dnapars.

        :param outfile:  Standard output filename.
        :param outtree:  Name of maximum parsimony tree.
        """
        infile = self._temp_infile(infile=self.infile)
        try:
            dnapars = pexpect.spawnu("dnapars %s" % infile)
            dnapars.sendline("Y\r")
            dnapars.waitnoecho()
            # TODO: Figure out how to catch output.
        except pexpect.EOF as e:
            self.phylip_log.exception(e)
        else:
            self._rename("outfile", outfile)
            self._rename("outtree", outtree)
        finally:
            os.remove(infile)

    def dnaml(self, outfile, outtree):
        """Generate a maximum likelihoood tree using dnapaml.

        :param outfile:  Standard output filename.
        :param outtree:  Name of maximum likelihoood tree.
        """
        infile = self._temp_infile(infile=self.infile)
        try:
            dnaml = pexpect.spawnu("dnaml %s" % infile)
            dnaml.sendline("Y\r")
            dnaml.waitnoecho()
        except pexpect.EOF as e:
            self.phylip_log.exception(e)
        else:
            self._rename("outfile", outfile)
            self._rename("outtree", outtree)
        finally:
            os.remove(infile)

    def dnadist(self, outfile):
        """Generate a distance matrix using dnadist.

        :param outfile:  distance matrix output filename.
        """
        infile = self._temp_infile(infile=self.infile)
        try:
            dnadist = pexpect.spawnu("dnadist %s" % infile)
            dnadist.sendline("Y\r")
            dnadist.waitnoecho()
        except pexpect.EOF as e:
            self.phylip_log.exception(e)
        else:
            self._rename("outfile", outfile)
        finally:
            os.remove(infile)
