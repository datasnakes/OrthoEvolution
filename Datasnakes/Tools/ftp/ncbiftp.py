"""Class to access NCBI's ftp server and easily download databases."""
import tarfile
import re
from time import time
from datetime import datetime
from multiprocessing.pool import ThreadPool
import os

from Datasnakes.Tools.ftp.baseftp import BaseFTPClient


class NcbiFTPClient(BaseFTPClient):
    """Access NCBI's FTP servers with ease."""
    def __init__(self, email):
        _NCBI = 'ftp.ncbi.nlm.nih.gov'
        super().__init__(_NCBI, email, keepalive=False, debug_lvl=0)
        self.blastdb_path = '/blast/db/'
        self.blastfasta_path = '/blast/db/FASTA/'
        self.refseqrelease_path = '/refseq/release/'
        self._taxdb = 'taxdb.tar.gz'  # Located in self.blastdb_path

        # TODO Use python to get these and turn into a json file or dict
        self.blastdbs = []
        self.blastfastadbs = []

    @classmethod
    def _pathformat(cls, path):
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

    @classmethod
    def extract_file(cls, file2extract):
        """Extract a tar.gz file."""
        if str(file2extract).endswith('tar.gz'):
            tar = tarfile.open(file2extract)
            tar.extractall()
            tar.close()
            print('Files were successfully extracted from %s.' % file2extract)
        else:
            raise ValueError('%s does not end in tar.gz' % file2extract)

    def listfiles(self, path='/'):
        """List all files in a path."""
        self._pathformat(path)
        _, files = self._walk(path)
        return files

    def listdirectories(self, path='/'):
        """List all directories in a path."""
        self._pathformat(path)
        directories, _ = self._walk(path)
        return directories

    def getrefseqrelease(self, database_name, data_type, file_type, download_path='', 
                         extract=True):
        """Download the preformatted blast database."""
        self.ftp.cwd(self.refseqrelease_path)
        releasedirs = self.listdirectories(self.refseqrelease_path)
        
        # Change to directory input
        if database_name not in releasedirs:
            raise FileNotFoundError('%s does not exist.' % database_name)
        
        self.ftp.cwd(database_name)
        releasefiles = self.listfiles()

        files2download = []
        for releasefile in releasefiles:
            if data_type and file_type in str(releasefile):
                files2download.append(releasefile)

        print('You are about to download theses files: %s' % files2download)

        # Download the files using multiprocessing
        download_time_secs = time()
        with ThreadPool(1) as download_pool:
            download_pool.map(self.download_file, files2download)
            minutes = (time() - download_time_secs) / 60
        print("Took %s minutes to download the files." % minutes)

        if extract:
            extract_time_secs = time()
            with ThreadPool(1) as extract_pool:
                extract_pool.map(self._extract_file, files2download)
                minutes = (time() - extract_time_secs) / 60
            print("Took %s minutes to extract from all files." % minutes)

    def getblastfasta(self, database_name, extract=True):
        """Download the fasta sequence database (not formatted)."""
        if str(database_name).startswith('est'):
            raise NotImplementedError('Est dbs cannot be downloaded yet.')
        self.ftp.cwd(self.blastfasta_path)
        blastfastafiles = self.listfiles(self.blastfasta_path)

        files2download = []
        for dbfile in blastfastafiles:
            if database_name in str(dbfile):
                files2download.append(dbfile)

        # Append the taxonomy database
        files2download.append(self._taxdb)
        print('You are about to download theses files: %s' % files2download)

        # Download the files using multiprocessing
        download_time_secs = time()
        with ThreadPool(1) as download_pool:
            download_pool.map(self.download_file, files2download)
            minutes = (time() - download_time_secs) / 60
        print("Took %s minutes to download the files." % minutes)

        if extract:
            extract_time_secs = time()
            with ThreadPool(1) as extract_pool:
                extract_pool.map(self._extract_file, files2download)
                minutes = (time() - extract_time_secs) / 60
            print("Took %s minutes to extract from all files." % minutes)

    def getblastdb(self, database_name, extract=True):
        """Download the fasta sequence database (not formatted)."""
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
        download_time_secs = time()
        with ThreadPool(1) as download_pool:
            download_pool.map(self.download_file, files2download)
            minutes = (time() - download_time_secs) / 60
        print("Took %s minutes to download the files." % minutes)

        if extract:
            extract_time_secs = time()
            with ThreadPool(1) as extract_pool:
                extract_pool.map(self._extract_file, files2download)
                minutes = (time() - extract_time_secs) / 60
            print("Took %s minutes to extract from all files." % minutes)

    def updatedb(self, database_path=os.getcwd(), update_days=7):
        """Check for when the database was last updated."""
        # Get a list of the files in the path
        filesinpath = os.listdir(database_path)
        for fileinpath in filesinpath:
            if str(fileinpath).endswith('.nal'):
                nalfile = str(fileinpath)
                dbname, ext = nalfile.split('.')
                del ext
                filetime = datetime.fromtimestamp(os.path.getctime(nalfile))
                format_filetime = filetime.strftime("%b %d, %Y at %I:%M:%S %p")

        print("Your database was last updated on: %s" % format_filetime)

        time_elapsed = datetime.now() - filetime
        if time_elapsed.days >= update_days:
            print('\nYour database needs updating.')
            self.getblastdb(dbname, download_path=database_path, extract=True)

        else:
            print('\nYour database has been updated within the last week.')
