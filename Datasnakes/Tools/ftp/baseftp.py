"""Base FTP module for connecting to FTP servers."""
from ftplib import FTP
import os

from Datasnakes.Tools.utils import FunctionRepeater


class BaseFTPClient(object):
    """The BaseFTP class provides basic functions for managing ftp clients."""
    def __init__(self, ftpsite, email, keepalive=False, debug_lvl=0):
        """Sets up the class by using the ftp site and email."""
        self._ftpsite = ftpsite
        self._email = email
        self.debug_lvl = debug_lvl
        self.ftp = self._login()
        self.keepalive = keepalive

        if self.keepalive:
            self._voidcmd_repeat, self._filetransfer_repeat = self._keepalive()

    def _keepalive(self):
        """Check to see if the FTP connection still exists using ftp.voidcmd

        Also, creates a filetransfer to prevent file transfer timeout.
        """
        voidcmd = FunctionRepeater(5, self.ftp.voidcmd, 'NOOP')
        filetransfer = FunctionRepeater(5, self._filetransfer, 'README.ftp')

        return voidcmd, filetransfer

    def _login(self):
        """Connect to the FTP server anonymously."""
        ftp = FTP(self._ftpsite, timeout=None)
        ftp.login(user='anonymous', passwd=self._email)
        ftp.set_debuglevel(self.debug_lvl)
        return ftp

    def close_connection(self):
        """Close the ftp connection."""
        self.ftp.close()

        if self.keepalive:
            self._voidcmd_repeat.stop()
            self._filetransfer_repeat.stop()

    def _filetransfer(self, filename):
        """Mimics a file transfer to prevent transfer time out."""
        print("WAIT: Keeping connection open by transferring file.")
        current_path = self.ftp.pwd()
        self.ftp.cwd('/')

        with open(filename, 'wb') as local_file:
            self.ftp.retrbinary('RETR %s' % filename, local_file.write)

        os.remove(filename)
        self.ftp.cwd(current_path)  # Return to starting path
        print("You may now resume your work.")
