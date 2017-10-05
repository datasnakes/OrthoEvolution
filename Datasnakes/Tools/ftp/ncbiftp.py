"""Class to access NCBI's ftp server and easily download databases."""
from baseftp import BaseFTPClient
import tarfile

from Datasnakes.Tools.parallel import Multiprocess


class NcbiFTPClient(BaseFTPClient):
    """Access NCBI's FTP servers with ease."""
    def __init__(self, email):
        _NCBI = 'ftp.ncbi.nlm.nih.gov'
        super().__init__(ftpsite=_NCBI, email=email, keepalive=False, debug_lvl=0)
        self.blastdb_path = '/blast/db/'
        self.blastfasta_path = '/blast/db/FASTA'
        blastdbs = []
        blastfastadbs = []

    class get:
        def blastdb(self, database_name, extract=True):
            """Download the preformatted blast database."""
            pass

        def blastfasta(self, database_name, extract=True):
            """Download the fasta sequence database (not formatted)."""

    def _extract(self):
        """Unzip blast databases."""
        pass

    def _download_file(self, filename):
        """Download the files one by one."""
        self.ftp.retrbinary("RETR " + filename, open(filename, 'wb').write)

    def _multi_download(self):
        """Download the files using multiprocessing."""
        Multiprocess(num_procs=4, function=self._download_file,
                     listinput='')

    def _listall(self, path):
        if path == 'blastdb':
            self.ftp.cwd(self.blastdb_path)
            blastdbs
            return
        elif path == 'blastfasta':
            self.ftp.cwd(self.blastdb_path)
        else:
            raise Exception

        filenames = self.ftp.nlst()

    def _parseMLSD(self):
        pass
