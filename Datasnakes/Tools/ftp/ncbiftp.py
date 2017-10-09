"""Class to access NCBI's ftp server and easily download databases."""
import tarfile
import re
import concurrent.futures
from time import time

from baseftp import BaseFTPClient
from Datasnakes.Tools.parallel import num_procs


class NcbiFTPClient(BaseFTPClient):
    """Access NCBI's FTP servers with ease."""
    def __init__(self, email):
        _NCBI = 'ftp.ncbi.nlm.nih.gov'
        super().__init__(ftpsite=_NCBI, email=email, keepalive=False, debug_lvl=0)
        self.blastdb_path = '/blast/db/'
        self.blastfasta_path = '/blast/db/FASTA/'
        self._taxdb = 'taxdb.tar.gz'  # Located in self.blastdb_path
        self._cpus = num_procs

        # Use python to get these and turn into a json file or dict
        self.blastdbs = []
        self.blastfastadbs = []

    def _pathformat(self, path):
        """Ensure proper formatting of the path."""
        pattern = re.compile('^/(.*?)/$')
        if not re.match(pattern, path):
            raise ValueError('Your path is not in a proper format.')

    def _walk(self, path):
        """Walk the ftp server and get files and directories."""
        file_list, dirs, nondirs = [], [], []
        try:
            self.ftp.cwd(path)
        except Exception as exp:
            print("Current path: ", self.ftp.pwd(), exp.__str__(), path)
            return [], []
        else:
            self.ftp.retrlines('LIST', lambda x: file_list.append(x.split()))
            for info in file_list:
                ls_type, name = info[0], info[-1]
                if ls_type.startswith('d'):
                    dirs.append(name)
                else:
                    nondirs.append(name)
            return dirs, nondirs

    def download_file(self, filename):
        """Download the files one by one."""
        with open(filename, 'wb') as localfile:
            self.ftp.retrbinary('RETR %s' % filename, localfile.write)
            print('%s was downloaded.' % str(filename))

    def extract_file(self, file2extract):
        if str(file2extract).endswith('tar.gz'):
            tar = tarfile.open(file2extract)
            tar.extractall()
            tar.close()
            print('')
        else:
            raise ValueError('%s does not end in tar.gz' % file2extract)

    def listfiles(self, path='/'):
        """List all files in a path."""
        self._pathformat(path)
        directories, files = self._walk(path)
        return files

    def listdirectories(self, path='/'):
        """List all directories in a path."""
        self._pathformat(path)
        directories, files = self._walk(path)
        return directories

    def getblastdb(self, database_name, download_path='', extract=False):
        """Download the preformatted blast database."""
        if str(database_name).startswith('est'):
            raise NotImplementedError('Est dbs cannot be downloaded yet.')
        self.ftp.cwd(self.blastdb_path)
        blastdbfiles = self.listfiles(self.blastdb_path)

        files2download = []
        for dbfile in blastdbfiles:
            if database_name in str(dbfile):
                files2download.append(dbfile)

        # Append the taxonomy database
        files2download.append(self._taxdb)
        print('You are about to download theses files: %s' % files2download)

        # Download the files using multiprocessing
        procs = self._cpus
        download_time_secs = time()
        with concurrent.futures.ProcessPoolExecutor(max_workers=procs) as executor:
            executor.map(self.download_file, files2download)
            minutes = (time() - download_time_secs) / 60
        print("Took %s minutes to download the files." % minutes)

        if extract:
            extract_time_secs = time()
            with concurrent.futures.ProcessPoolExecutor(max_workers=procs) as executor:
                executor.map(self._extract_file, files2download)
                minutes = (time() - extract_time_secs) / 60
            print("Took %s minutes to extract from all files." % minutes)

    def getblastfasta(self, database_name, extract=True):
        """Download the fasta sequence database (not formatted)."""
        self.ftp.cwd(self.blastdb_path)
        blastdbfiles = self.listfiles(self.blastdb_path)

        files2download = []
        for dbfile in blastdbfiles:
            if database_name in str(dbfile):
                files2download.append(dbfile)

        # Append the taxonomy database
        files2download.append(self._taxdb)
        print('You are about to download theses files: %s' % files2download)

        downloaded = []
        for file in files2download:
            self.download_file(file)
            downloaded.append(file)

        if extract:
            for downloaded_file in downloaded:
                self.extract_file(downloaded_file)
