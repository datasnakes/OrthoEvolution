"""Base FTP client for connecting to remote file repositories."""

import logging
from ftplib import FTP, all_errors
from pathlib import Path, PurePosixPath
from tempfile import TemporaryFile

from OrthoEvol.utilities import FunctionRepeater


logger = logging.getLogger(__name__)


class BaseFTPClient:
    """Provide shared connection management for FTP clients."""

    def __init__(
        self,
        ftpsite: str,
        user: str,
        password: str,
        keepalive: bool = False,
        debug_lvl: int = 0,
        timeout: float = 600.0,
    ) -> None:
        """Connect to an FTP server with the supplied credentials.

        :param ftpsite: Host name of the FTP server.
        :param user: User name used to log in.
        :param password: Password used to log in.
        :param keepalive: Whether to periodically keep the connection active.
        :param debug_lvl: Verbosity level for the FTP connection.
        :param timeout: Socket timeout in seconds.
        """
        self._ftpsite = ftpsite
        self._user = user
        self._password = password
        self._debug_lvl = debug_lvl
        self._timeout = timeout
        self.ftp = self._login()
        self.__keepalive = keepalive

        if self.__keepalive:
            self._voidcmd_repeat, self._filetransfer_repeat = self._keepalive()

    def _keepalive(self) -> tuple[FunctionRepeater, FunctionRepeater]:
        """Create periodic commands that keep an idle FTP session active.

        .. warning:: :func:`_keepalive` is not well tested.
           Avoid using it if possible.
        """
        voidcmd = FunctionRepeater(5, self.ftp.voidcmd, "NOOP")
        filetransfer = FunctionRepeater(5, self._filetransfer, "README.ftp")
        return voidcmd, filetransfer

    def _login(self) -> FTP:
        """Open and validate an FTP connection."""
        ftp = FTP(self._ftpsite, timeout=self._timeout)
        ftp.login(user=self._user, passwd=self._password)
        ftp.voidcmd("NOOP")
        ftp.set_debuglevel(self._debug_lvl)
        return ftp

    def close_connection(self) -> None:
        """Stop keepalive workers and close the FTP connection."""
        if self.__keepalive:
            self._voidcmd_repeat.stop()
            self._filetransfer_repeat.stop()

        try:
            self.ftp.quit()
        except all_errors:
            # A timed-out server may reject QUIT but still hold a local socket.
            self.ftp.close()

    def _filetransfer(self, filename: str | Path) -> None:
        """Transfer and discard a small file to keep the session active.

        :param filename: Remote path of the file used for the keepalive transfer.
        """
        remote_path = PurePosixPath(str(filename))
        current_path = self.ftp.pwd()
        logger.info("Keeping the FTP connection active with %s.", remote_path)

        try:
            self.ftp.cwd("/")
            with TemporaryFile() as temporary_file:
                self.ftp.retrbinary(
                    f"RETR {remote_path.as_posix()}",
                    temporary_file.write,
                )
        finally:
            self.ftp.cwd(current_path)
