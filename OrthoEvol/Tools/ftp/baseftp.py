"""Base FTP module for connecting to FTP servers."""
from ftplib import FTP, error_perm
import os
import contextlib

from OrthoEvol.utilities import FunctionRepeater


class BaseFTPClient(object):
    """The BaseFTP class provides basic functions for managing ftp clients.

    .. seealso:: :class:`NcbiFTPClient`
    """

    def __init__(self, ftpsite, user, password, keepalive=False, debug_lvl=0):
        """Connect to a ftp site using a username and password.

        :param ftpsite: The url or http address of the ftp site you want to connect to.
        :param user: The name of the user that will log in.
        :type user: str
        :param password: The password needed to log in to the ftp site.
        :type password: str
        :param keepalive: Flag to determine whether to keepalive the connection, defaults to False
        :type keepalive: bool, optional
        :param debug_lvl: Verbosity level for debugging ftp connection, defaults to 0
        :type debug_lvl: int, optional
        """
        self._ftpsite = ftpsite
        self._user = user
        self._password = password
        self._debug_lvl = debug_lvl
        self.ftp = self._login()
        self.__keepalive = keepalive

        if self.__keepalive:
            self._voidcmd_repeat, self._filetransfer_repeat = self._keepalive()

    def _keepalive(self):
        """Check to see if the FTP connection still exists using ftp.voidcmd

        Also, creates a filetransfer to prevent file transfer timeout.

        .. warning:: :func:`_keepalive` is not well tested.
        Avoid using it if possible.
        """

        voidcmd = FunctionRepeater(5, self.ftp.voidcmd, 'NOOP')
        filetransfer = FunctionRepeater(5, self._filetransfer, 'README.ftp')

        return voidcmd, filetransfer

    def _login(self):
        """Connect to the FTP server."""

        with contextlib.suppress(error_perm):
            ftp = FTP(self._ftpsite, timeout=600)
            ftp.login(user=self._user, passwd=self._password)
            ftp.voidcmd('NOOP')
            ftp.set_debuglevel(self._debug_lvl)

        return ftp

    def close_connection(self):
        """Close the ftp connection."""

        self.ftp.close()

        if self.__keepalive:
            self._voidcmd_repeat.stop()
            self._filetransfer_repeat.stop()

    def _filetransfer(self, filename):
        """Mimics a file transfer to prevent transfer time out.

        :param filename: Path of file to download.
        """

        print("WAIT: Keeping connection open by transferring file.")
        current_path = self.ftp.pwd()
        self.ftp.cwd('/')

        with open(filename, 'wb') as local_file:
            self.ftp.retrbinary('RETR %s' % filename, local_file.write)

        os.remove(filename)
        self.ftp.cwd(current_path)  # Return to starting path
        print("You may now resume your work.")
