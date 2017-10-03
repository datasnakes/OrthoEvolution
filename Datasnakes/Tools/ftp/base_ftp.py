"""FTP module for"""
from ftplib import FTP
# import time
# import shutil
# import subprocess
# import configparser


class BaseFTP(object):
    """The BaseFTP class provides basic functions for managing ftp clients."""
    def __init__(self, ftpsite, email):
        """Initialize the class by using the ftp site and email."""
        self.site = ftpsite
        self.email = email

    def ftp_check(self):
        # TODO reevaluate how to do this.
        """Check to see if the FTP connection still exists.

        If it doesn't then it reconnects and returns an instance of the
        connection.
        """
        ftp = self.ftp_connect(self.__NCBI_FTP, self.email)
        ftp.voidcmd('NOOP')
        ftp.cwd(self.path)
        return ftp

    @classmethod
    def ftp_login(self):
        """Connect to the FTP server anonymously."""
        ftp = FTP(self.site, timeout=None)
        ftp.login(user='anonymous', passwd=self.email)
        return ftp
