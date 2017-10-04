"""Class to access NCBI's ftp server and easily download databases."""
from baseftp import BaseFTPClient


class NcbiFTPClient(BaseFTPClient):
    """Access NCBI's FTP servers with ease."""
    def __init__(self, email):
        _NCBI = 'ftp.ncbi.nlm.nih.gov'
        super().__init__(ftpsite=_NCBI, email=email, keepalive=False, debug_lvl=0)

    def list_blastdbs(self):
        # TODO Create function to list databases
        """List all blast databases."""
        pass

    def getblastdb(self, database_name, unzip=True):
        # TODO Create function to download blast databases
        """Download the blast database."""
        pass

    def _unzip(self):
        """Unzip blast databases."""
        pass
